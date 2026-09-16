
assign_letters <- function(contrasts_df, data, y_var, group_var, offset = 0.3, alpha = 0.05) {
  groups <- unique(c(
    sub(" - .*", "", contrasts_df$contrast),
    sub(".* - ", "", contrasts_df$contrast)
  ))
  
  ns_pairs <- contrasts_df %>%
    filter(p.value > alpha) %>%
    mutate(g1 = sub(" - .*", "", contrast),
           g2 = sub(".* - ", "", contrast))
  
  letter_df <- data.frame(group = groups, letter = NA)
  current_letter <- 1
  
  for (g in groups) {
    shared <- ns_pairs$g2[ns_pairs$g1 == g]
    existing <- letter_df$letter[letter_df$group %in% shared]
    existing <- existing[!is.na(existing)]
    
    if (length(existing) > 0) {
      letter_df$letter[letter_df$group == g] <- existing[1]
    } else {
      letter_df$letter[letter_df$group == g] <- letters[current_letter]
      current_letter <- current_letter + 1
    }
  }
  
  letter_df[[group_var]] <- as.numeric(gsub("[^0-9.]", "", letter_df$group))
  letter_df <- letter_df[, c(group_var, "letter")]  # drop group column before join
  
  y_positions <- data %>%
    group_by(across(all_of(group_var))) %>%
    summarise(y_pos = log(max(.data[[y_var]], na.rm = TRUE)) + offset, .groups = "drop")
  
  letter_df <- left_join(letter_df, y_positions, by = group_var)
  return(letter_df)
}
