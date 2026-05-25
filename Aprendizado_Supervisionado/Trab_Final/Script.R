library(tidyverse)
library(tidytext)

# Poemas de Shakespeare:
poemas_shake <- Rpoet::get.poetry("author", "William Shakespeare", "lines")

poemas_chaucer <- Rpoet::get.poetry("author", "Geoffrey Chaucer", "lines")

# Léxico:
bing <- get_sentiments("bing")

sentimentos_shake <- vector("integer", dim(poemas_shake)[1])

for (i in 1:dim(poemas_shake)[1]) {
    palavras_shake <- tibble(text = poemas_shake[[1]][[i]]) %>%
        unnest_tokens(output = "word", input = "text") %>%
        anti_join(get_stopwords())
    sentimentos_shake[i] <- palavras_shake %>%
        inner_join(bing, relationship = "many-to-many") %>%
        rownames_to_column(var = "linha") %>%
        mutate(sentiment = case_when(sentiment == "negative" ~ -1,
                                     .default = 1)) %>%
        count(wt = sentiment) %>%
        unlist()
}

sentimentos_chaucer <- vector("integer", dim(poemas_chaucer)[1])

for (i in 1:dim(poemas_chaucer)[1]) {
    palavras_chaucer <- tibble(text = poemas_chaucer[[1]][[i]]) %>%
        unnest_tokens(output = "word", input = "text") %>%
        anti_join(get_stopwords())
    sentimentos_chaucer[i] <- palavras_chaucer %>%
        inner_join(bing, relationship = "many-to-many") %>%
        rownames_to_column(var = "linha") %>%
        mutate(sentiment = case_when(sentiment == "negative" ~ -1,
                                     .default = 1)) %>%
        count(wt = sentiment) %>%
        unlist()
}



tibble(sent_shake = sentimentos_shake) %>%
    arrange(sent_shake) %>%
    mutate(linha = row_number()) %>%
    ggplot(aes(x = linha, y = sent_shake)) +
    geom_bar(stat = "identity")

tibble(sent_chaucer = sentimentos_chaucer) %>%
    arrange(sent_chaucer) %>%
    mutate(linha = row_number()) %>%
    ggplot(aes(x = linha, y = sent_chaucer)) +
    geom_bar(stat = "identity")
