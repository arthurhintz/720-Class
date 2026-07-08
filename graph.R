library(tidyverse)

dat <- read.table("results.txt")

dat <- dat |> 
  mutate(trt_eff = case_when(trt_id < 2 ~ 0,
                             trt_id < 6 ~ 1,
                             trt_id < 10 ~ 2,
                             trt_id < 13 ~ 3,
                             trt_id >= 13 ~ 4))


colnames(dat)


plot(dat$entropy, dat$field)


plot(dat$entropy, dat$block_n)
plot(dat$entropy, dat$trt_id)
plot(dat$upper_ci, dat$lower_ci)

dat
resu1 <- dat |> 
  group_by(trt_eff, block_n) |> 
  summarise(entropy = mean(entropy)) |> 
  pivot_wider(names_from = block_n, values_from = entropy, names_prefix = "block_") |>
  pivot_longer(
    cols = -c(trt_eff, block_4),
    names_to = "block_n",
    values_to = "entropy_outro_bloco"
  ) |>
  mutate(
    block_n = readr::parse_number(block_n),
    diff_entropy = block_4 - entropy_outro_bloco  # bloco 4 menos o outro bloco
  ) |>
  select(trt_eff, block_n, entropy_block4 = block_4, entropy_outro_bloco, diff_entropy) |>
  arrange(trt_eff, block_n)

resu1

resu1 |> 
  ggplot(aes(trt_eff, diff_entropy, colour = factor(block_n))) +
  geom_line() +
  geom_hline(yintercept = 0)



#==========/==========/==========/==========/==========/==========/==========/==========/

check_coverage <- function(media, ci_low, ci_up){
  as.integer(ci_low <= media & ci_up >= media)
}

resu2 <- dat |>
  filter(trt_eff != 0) |> 
  select(field, block_n, lower_ci, upper_ci, trt_eff)

resu2 |> count(trt_eff)

coverage_summary <- resu2 |>
  mutate(
    count = case_when(
      trt_eff == 1 ~ check_coverage(-2, lower_ci, upper_ci),
      trt_eff == 2 ~ check_coverage( 0, lower_ci, upper_ci),
      trt_eff == 3 ~ check_coverage( 2, lower_ci, upper_ci),
      trt_eff == 4 ~ check_coverage( 4, lower_ci, upper_ci)
    )
  ) |>
  group_by(block_n) |>
  summarise(
    n_total = n(),
    n_covered = sum(count),
    pct_covered = 100 * mean(count),
    .groups = "drop"
  ) |>
  mutate(Block = factor(block_n))

p2 <- ggplot(coverage_summary,
             aes(x = Block, y = pct_covered)) +
  geom_col(fill = "grey70", colour = "black", width = 0.7) +
  geom_hline(yintercept = 95, linetype = "dashed",
             colour = "red", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f", pct_covered)),
            vjust = -0.4, size = 3.5) +
  scale_y_continuous(
    limits = c(0, 105), breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Block",y = "Coverage (%)") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank())

p2


ggsave("figures/coverage.pdf", p2, dpi = 250, width = 9, height = 6)

#==========/==========/==========/==========/==========/==========/==========/==========/

resu2 <- dat |> 
  filter(trt_eff == 2) |> 
  mutate(Block = block_n) |>
  select(field,Block,mean_post)

media_b <- resu2 |> 
  group_by(Block) |> 
  summarise(media_b  = mean(mean_post),
            mse = mean(mean_post^2))


p3 <- ggplot(resu2, aes(x = mean_post)) +
  geom_density(alpha = 0.4, linewidth = 0.6) + 
  geom_vline(xintercept = 0, color = "red", linewidth = 0.5) + 
  geom_vline(
    data = media_b, aes(xintercept = media_b),
    linetype = "dashed", color = "black", linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = media_b,
    aes(x = media_b, y = Inf,
        label = sprintf("mean = %.2f\nMSE = %.2f", media_b, mse)),
    vjust = 1.3, hjust = -0.8, size = 3.2, lineheight = 0.9,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Block, ncol = 3, labeller = label_both) + 
  labs(
    x = expression("Posterior mean of " ~ T[list(6, "...", 9)]),
    y = "Density"
  ) +
  theme_bw(base_size = 12)
p3
ggsave("figures/density_trt.pdf", p3, dpi = 250, width = 9, height = 6)

#==========/==========/==========/==========/==========/==========/==========/==========/


resu2 <- dat |> 
  filter(trt_eff != 0) |> 
  mutate(
    Block = block_n,
    ci_width = upper_ci - lower_ci
  ) |> 
  select(field, Block, ci_width)

media_width <- resu2 |> 
  group_by(Block) |> 
  summarise(media_width = mean(ci_width), .groups = "drop")

p4 <- ggplot(resu2, aes(x = ci_width)) +
  geom_density(alpha = 0.4, linewidth = 0.6) + 
  geom_vline(
    data = media_width, aes(xintercept = media_width),
    linetype = "dashed", color = "black", linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = media_width,
    aes(x = media_width, y = Inf, label = sprintf("%.2f", media_width)),
    vjust = 1.5, hjust = -0.5, size = 3.2, inherit.aes = FALSE
  ) +
  facet_wrap(~ Block, ncol = 3, labeller = label_both) + 
  labs(
    x = "Credible interval width",
    y = "Density"
  ) +
  theme_bw(base_size = 12)
p4

ggsave("figures/CI_width.pdf", p4, dpi = 250, width = 9, height = 6)






#==========/==========/==========/==========/==========/==========/==========/==========/


resu3 <- dat |> 
  group_by(trt_id, block_n) |> 
  summarise(up_ic = mean(upper_ci),
            low_ic = mean(lower_ci),
            amp = up_ic - low_ic)

resu3


#==========/==========/==========/==========/==========/==========/==========/==========/

resu4 <- dat |>
  select(field, block_n, trt_id, entropy) |>
  pivot_wider(names_from = block_n, values_from = entropy, names_prefix = "block_") |>
  rowwise() |>
  mutate(
    entropy_block4 = block_4,
    min_outros = min(c_across(-c(field, trt_id, block_4)), na.rm = TRUE),
    counter = if_else(entropy_block4 <= min_outros, 1L, 0L)
  ) |>
  ungroup() |>
  select(field, trt_id, entropy_block4, min_outros, counter)

resu4

mean(resu4$counter)

resu4 |>
  group_by(trt_id) |>
  summarise(
    n_total  = n(),                
    n_wins   = sum(counter),
    entropy = mean(entropy_block4),
    pct_wins = 100 * mean(counter),
    .groups = "drop"
  )


# Em nenhum tratamento teve pelo um completed block desing que se destacou 
# em relacao a pelo menos 1 bloco do inclompleted desing



#==========/==========/==========/==========/==========/==========/==========/==========/


dat |> 
  ggplot(aes(y = entropy, group = factor(block_n), fill = factor(block_n))) +
  geom_boxplot()

#==========/==========/==========/==========/==========/==========/==========/==========/
#Power

resu5 <- dat |> 
  filter(trt_id != 1) |> 
  group_by(trt_eff, block_n) |> 
  summarise(power = round(mean(posterior_prob), 3))



resu5 |> 
  ggplot(aes(block_n, power, colour = factor(trt_eff))) +
  geom_point()




resu5 |> 
  ggplot(aes(trt_eff, power, colour = factor(block_n))) +
  geom_point()






