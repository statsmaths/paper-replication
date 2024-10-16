suppressPackageStartupMessages({
library(tidyverse)
library(stringdist)
library(stringi)
library(cleanNLP)
library(xtable)
})

options(dplyr.summarise.inform = FALSE)
options(width = 120)
theme_set(theme_minimal())
cnlp_init_udpipe()

read_whisper <- function(path)
{
  z <- read_delim(
    path,
    delim = "\t",
    quote = "",
    col_types = "cddd",
    col_names = c("token", "score", "xmin", "xmax")
  )
  z$xmin <- z$xmin * 10
  z$xmax <- z$xmax * 10
  z  
}

# create a version of the whisper data that removes special tokens
# and has data about the entire subtoken; this is the dataset that
# we will be matching the reference text into
talign_api_clean <- function(api)
{
  api_subtokens <- filter(api, stri_sub(token, 1L, 1L) != "[") |>
    mutate(is_start = !(stri_sub(token, 1L, 1L) %in% c(LETTERS, letters))) |>
    mutate(token_id = cumsum(is_start)) |>
    group_by(token_id) |>
    mutate(
      n_subtoken = n(),
      full = stri_trim(paste(token, collapse = "")),
      full_alt = stri_trim(paste(token, collapse = "|"))
    ) |>
    group_by()  

  return(api_subtokens)
}

# create a data frame has one row per full token from the API
talign_api_token <- function(api_subtokens)
{
  api_token <- api_subtokens |>
    group_by(token_id) |>
    summarize(
      full = first(full),
      full_alt = first(full_alt),
      score = first(score),
      xmin = min(xmin),
      xmax = max(xmax)
    )

  return(api_token)
}

# produce a tokenized version of the reference text; need to
# "fix" the way the quotes are done to match the way that
# OpenAI does tokenization
talign_ref_tokens <- function(ref)
{
  token_ref <- cnlp_annotate(ref)$token$token
  while(any(token_ref == "\""))
  {
    idx <- which(token_ref == "\"")
    token_ref[idx + 1L] <- sprintf("\"%s", token_ref[idx + 1L])
    token_ref <- token_ref[-idx]
  }
  return(token_ref)
}

# here is the core of the algorithm; cycle through the elements
# of v1 looking for matches in v2; if we go to far without a 
# match consider this a missing token and go to the next; we
# keep track of everything relative to the index of each sequence
# through the vectors k1 and k2.
talign_alignment <- function(v1, v2, md = 3L)
{
  # find unique terms that do match; these are anchors for the 
  # text to avoid making large errors
  wunique <- setdiff(v1, v1[duplicated(v1)])
  wunique <- setdiff(v1, v2[duplicated(v2)])
  wunique <- wunique[wunique %in% v2]

  idx1 <- match(wunique, v1)
  idx2 <- match(wunique, v2)

  del1 <- diff(idx1)
  del2 <- diff(idx2)
  index <- which(abs(del1 - del2) > 5L | del2 < 0)
  if (length(index) > 0)
  {
    idx1 <- idx1[-index]
    idx2 <- idx2[-index]
  }

  # v_index is a vector following the reference text with pointers
  # into the transcription pointing to the best alignment; k1 is a
  # scalar pointer to the current element of the reference text and
  # k2 is a scalar pointer to the current element in the
  # transcription text
  v_index <- rep(NA_integer_, length(v1))
  v_index[idx1] <- idx2

  k1 <- 1L
  k2 <- 1L
  while(TRUE)
  {
    if (!is.na(v_index[k1]))
    {
      k2 <- v_index[k1] + 1L
      k1 <- k1 + 1L
    } else {
      if (v1[k1] == v2[k2])
      {
        v_index[k1] <- k2
        k1 <- k1 + 1L
        k2 <- k2 + 1L
      } else {
        mset <- (v1[k1] == v2[seq(k2 + 1, k2 + md)])
        if (any(mset, na.rm = TRUE))
        {
          k2 <- k2 + which.max(mset)
          v_index[k1] <- k2
          k1 <- k1 + 1L
          k2 <- k2 + 1L
        } else {
          k1 <- k1 + 1L
        }        
      }
    }
    if ((k1 > length(v1) | (k2 > length(v2)))) { break }
  }

  # now, in the output let's put this back together into a data frame
  # with one row for each token in the reference text
  df <- tibble(
    ref_text = v1,
    api_text = NA_character_,
    api_text_alt = NA_character_,
    score = NA_real_,
    xmin = NA_real_,
    xmax = NA_real_,
    id = v_index,
    block_id = cumsum(!is.na(if_else(!is.na(lag(v_index)), 1L, v_index)))
  ) |>
    group_by(block_id) |>
    mutate(gap_size_ref = n()) |>
    ungroup()

  # show where API text has gaps relative to the reference text
  df$api_text[!is.na(v_index)] <- v2[v_index[!is.na(v_index)]]
  df$api_text_alt[!is.na(v_index)] <- api_token$full_alt[v_index[!is.na(v_index)]]
  df$score[!is.na(v_index)] <- api_token$score[v_index[!is.na(v_index)]]
  df$xmin[!is.na(v_index)] <- api_token$xmin[v_index[!is.na(v_index)]]
  df$xmax[!is.na(v_index)] <- api_token$xmax[v_index[!is.na(v_index)]]

  # construct token ids to match back to the API data
  df$c1 <- cummax(if_else(is.na(df$id), -Inf, df$id))
  df$c2 <- rev(cummin(rev(if_else(is.na(df$id), Inf, df$id))))
  df$c1[is.na(df$api_text)] <- df$c1[is.na(df$api_text)] + 1L
  df$c2[is.na(df$api_text)] <- df$c2[is.na(df$api_text)] - 1L
  df$gap_size_api <- df$c2 - df$c1
  idx <- which((is.na(df$id)) & is.finite(df$gap_size_api))
  df$gap_text <- ""
  df$gap_text[idx] <- map2_chr(df$c1[idx], df$c2[idx], ~ paste(v2[seq(..1, ..2)], collapse = " "))

  # the above doesn't work when there is no corresponding API text;
  # fix these here
  df$c1[df$c1 > df$c2] <- NA_integer_
  df$c2[is.na(df$c1)] <- NA_integer_
  df$gap_size_api[is.na(df$c1)] <- 0
  df$gap_text[is.na(df$c1)] <- ""

  return(df)
}

# this function takes the raw alignement, which has exactly one row
# for each token in the reference text, and breaks it down into blocks
# of matching tokens; this is better for doing metrics on the comparison
talign_alignment_clean <- function(align_raw, api_token)
{
  idx <- which(is.na(align_raw$api_text_alt) & !duplicated(align_raw$block_id))
  align_raw$api_text[idx] <- stri_replace_all(align_raw$gap_text[idx], "", fixed = "|")
  align_raw$api_text_alt[idx] <- align_raw$gap_text[idx]
  align_raw$api_text[is.na(align_raw$api_text)] <- ""
  align_raw$api_text_alt[is.na(align_raw$api_text_alt)] <- ""

  idx <- which(is.na(align_raw$score))
  align_raw$score[idx] <- api_token$score[align_raw$c1[idx]]
  align_raw$xmin[idx] <- api_token$xmin[align_raw$c1[idx]]
  align_raw$xmax[idx] <- api_token$xmax[align_raw$c2[idx]]

  ppt <- function(x) { return(
    stri_replace_all(stri_trans_tolower(x), "", fixed = " ")
  ) } 

  align <- align_raw |>
    group_by(block_id) |>
    summarize(
      ref_text = paste(ref_text, collapse = " "),
      api_text = paste(api_text, collapse = " "),
      api_text_alt = paste(api_text_alt, collapse = " "),
      ref_size = n(),
      api_size = first(c2 - c1) + 1,
      api_size_sub = api_size + stri_count(api_text_alt, fixed = "|"),
      sdist = stringdist(ppt(ref_text), ppt(api_text)),
      score = first(score),
      xmin = min(xmin),
      xmax = max(xmax)
    )

  align$api_size[is.na(align$api_size)] <- 0
  align$api_size_sub[is.na(align$api_size_sub)] <- 0
  return(align)
}


talign_alignment2 <- function(v1, v2, md = 3L)
{
  # find unique terms that do match; these are anchors for the 
  # text to avoid making large errors
  wunique <- setdiff(v1, v1[duplicated(v1)])
  wunique <- setdiff(v1, v2[duplicated(v2)])
  wunique <- wunique[wunique %in% v2]

  idx1 <- match(wunique, v1)
  idx2 <- match(wunique, v2)

  del1 <- diff(idx1)
  del2 <- diff(idx2)
  index <- which(abs(del1 - del2) > 5L | del2 < 0)
  if (length(index) > 0)
  {
    idx1 <- idx1[-index]
    idx2 <- idx2[-index]
  }

  # v_index is a vector following the reference text with pointers
  # into the transcription pointing to the best alignment; k1 is a
  # scalar pointer to the current element of the reference text and
  # k2 is a scalar pointer to the current element in the
  # transcription text
  v_index <- rep(NA_integer_, length(v1))
  v_index[idx1] <- idx2

  k1 <- 1L
  k2 <- 1L
  while(TRUE)
  {
    if (!is.na(v_index[k1]))
    {
      k2 <- v_index[k1] + 1L
      k1 <- k1 + 1L
    } else {
      if (v1[k1] == v2[k2])
      {
        v_index[k1] <- k2
        k1 <- k1 + 1L
        k2 <- k2 + 1L
      } else {
        mset <- (v1[k1] == v2[seq(k2 + 1, k2 + md)])
        if (any(mset, na.rm = TRUE))
        {
          k2 <- k2 + which.max(mset)
          v_index[k1] <- k2
          k1 <- k1 + 1L
          k2 <- k2 + 1L
        } else {
          k1 <- k1 + 1L
        }        
      }
    }
    if ((k1 > length(v1) | (k2 > length(v2)))) { break }
  }

  # now, in the output let's put this back together into a data frame
  # with one row for each token in the reference text
  df <- tibble(
    ref_text = v1,
    api_text = NA_character_,
    api_text_alt = NA_character_,
    score = NA_real_,
    xmin = NA_real_,
    xmax = NA_real_,
    id = v_index,
    block_id = cumsum(!is.na(if_else(!is.na(lag(v_index)), 1L, v_index)))
  ) |>
    group_by(block_id) |>
    mutate(gap_size_ref = n()) |>
    ungroup()

  # show where API text has gaps relative to the reference text
  df$api_text[!is.na(v_index)] <- v2[v_index[!is.na(v_index)]]
  #df$api_text_alt[!is.na(v_index)] <- api_token$full_alt[v_index[!is.na(v_index)]]
  #df$score[!is.na(v_index)] <- api_token$score[v_index[!is.na(v_index)]]
  #df$xmin[!is.na(v_index)] <- api_token$xmin[v_index[!is.na(v_index)]]
  #df$xmax[!is.na(v_index)] <- api_token$xmax[v_index[!is.na(v_index)]]

  # construct token ids to match back to the API data
  df$c1 <- cummax(if_else(is.na(df$id), -Inf, df$id))
  df$c2 <- rev(cummin(rev(if_else(is.na(df$id), Inf, df$id))))
  df$c1[is.na(df$api_text)] <- df$c1[is.na(df$api_text)] + 1L
  df$c2[is.na(df$api_text)] <- df$c2[is.na(df$api_text)] - 1L
  df$gap_size_api <- df$c2 - df$c1
  idx <- which((is.na(df$id)) & is.finite(df$gap_size_api))
  df$gap_text <- ""
  df$gap_text[idx] <- map2_chr(df$c1[idx], df$c2[idx], ~ paste(v2[seq(..1, ..2)], collapse = " "))

  # the above doesn't work when there is no corresponding API text;
  # fix these here
  df$c1[df$c1 > df$c2] <- NA_integer_
  df$c2[is.na(df$c1)] <- NA_integer_
  df$gap_size_api[is.na(df$c1)] <- 0
  df$gap_text[is.na(df$c1)] <- ""

  return(df)
}