# url <- "https://github.com/Ngyein-Min-Zaw/AMH/raw/refs/heads/main/Data/UV%20(04.08.2006)%20-%20completed.xlsx"
# download.file(url,
#      destfile = "./Data/UV (04.08.2006) - completed.xlsx",
#      mode = "wb"
# )


library(rsm)
library(tidyverse)
library(extrafont)
library(readxl)
library(ggtext)

## Load Data

# url <- "https://github.com/Ngyein-Min-Zaw/AMH/raw/refs/heads/main/Data/(26.06.2026)%20absorbance%20data.xlsx"
# download.file(url,
#               destfile = "./Data/(26.06.2026) absorbance data.xlsx",
#               mode = "wb"
# )
library(readxl)
UV <- read_excel("./Data/downloaded_file.xlsx")

RSM <- read_excel("./Data/(26.06.2026) absorbance data.xlsx")
RSM2 <- RSM %>% mutate(Temperature = str_remove(Temp, "°C") %>% as.numeric())

mfont <- "Times New Roman"


codedrsm2 <- coded.data(
     RSM2,
     x1 ~ (Temperature - 50) / 10,
     x2 ~ (pH - 6) / 0.5
)

print(codedrsm2, code = F)
rsmmod <- rsm(avg ~ SO(x1, x2), codedrsm2)
summary(rsmmod)
canonical(rsmmod)

## CCD Design figure
for_rsm_vis <- left_join(RSM2, codedrsm2)
for_rsm_vis |> ggplot(aes(x1, x2)) +
     geom_point() +
     geom_text(
          aes(x1, x2 + 0.25,
               label = paste0("(", Temp, ", ", pH, ")")
          ),
          family = "Times New Roman",
          size = 4
     ) +
     coord_cartesian(xlim = c(-2, 2), ylim = c(-2, 2)) +
     theme_minimal() +
     labs(
          x = "x1 ~ (Temperature - 50)/10",
          y = "x2 ~ (pH - 6)/0.5"
     ) +
     theme(
          text = element_text(
               family = "Times New Roman",
               size = 12
          )
     )
ggsave("./Figs/ccd_design.png", width = 4.3, height = 4, dpi = 300)

### RSM 3D plot ###
png(
     file = "./Figs/3D_rsm.png", width = 5.5, height = 5.5, units = "in",
     res = 100, pointsize = 12
)
persp(rsmmod, ~ x1 + x2,
     at = canonical(rsmmod)$xs,
     contours = list(z = "bottom"), ticktype = "detailed",
     # nticks = 7,
     # shade = .2,
     # phi = 20,
     col = hcl.colors(20, "Purple-Yellow"),
     family = "Times New Roman",
     zlab = "\n\nAbsorbance",
     xlab = c("\nTemperature", "\npH"),
     cex.lab = 1,
     cex.axis = 1
)
dev.off()
###################
## Enzyme Activity with Four Starch at varying concentration at different duration##


UV_04_08_2006_ <- read_excel("./Data/UV (04.08.2006) - completed.xlsx",
     sheet = "Sheet1"
)
UV <- UV_04_08_2006_[1:42, 1:6]
UV <- UV %>% filter(!is.na(`0 min`))

ns <- names(UV)
ns[1] <- "Starch"
names(UV) <- ns

UV <- UV %>% filter(!str_detect(`0 min`, "min"))

starch <- c("Soluble Starch", "Topioca Starch", "Corn Starch", "Rice Starch")
starch <- c(rep(starch, each = 2), "Negative control")
starch <- rep(starch, times = 4)
UV[, 1] <- starch

lvl <- c(
     "Soluble Starch",
     "Corn Starch",
     "Rice Starch",
     "Topioca Starch",
     "Negative control"
)

UV2 <- UV %>%
     mutate(sub_conc = rep(c(1, 0.5, 1.5, 2), each = 9)) %>%
     pivot_longer(2:6,
          names_to = "Time",
          values_to = "ABS"
     ) %>%
     mutate(
          Time = str_remove(Time, " min"),
          Time = as.double(Time)
     ) %>%
     mutate(
          Starch = factor(Starch,
               levels = lvl
          ),
          ABS = as.numeric(ABS)
     )

## plot side by side ##
f_starch <- UV2 %>%
     # filter(Starch!="Negative control")%>%
     ggplot(aes(x = Time, y = ABS, color = Starch)) +
     geom_point() +
     coord_cartesian(ylim = c(0, 2)) +
     geom_smooth(
          method = "lm", se = F,
          formula = y ~ poly(x, 2, raw = T)
     ) +
     scale_x_continuous(breaks = seq(0, 60, 20)) +
     scale_y_continuous(breaks = seq(0, 2, 0.5)) +
     scale_color_manual(
          breaks = lvl,
          values = c(
               "green",
               "yellow",
               "black",
               "blue",
               "red"
          )
     ) +
     guides(color = guide_legend(nrow = 2, byrow = TRUE)) + # ← key line
     labs(
          x = "Time (minutes)",
          y = "Absorbance at 620 nm",
          color = ""
     ) +
     facet_grid(. ~ sub_conc,
          labeller = as_labeller(function(x) {
               paste(x, "%")
          })
     ) +
     theme_bw() +
     theme(
          text = element_text(
               family = mfont,
               size = 12,
               color = "black"
          ),
          legend.position = "bottom",
          axis.text = element_text(
               size = 12,
               color = "black"
          ),
          axis.title = element_text(
               size = 12,
               color = "black"
          ),
          axis.line = element_line(color = "black"),
          legend.text = element_text(size = 12, family = mfont),
          strip.background = element_rect(fill = "grey90"),
          strip.text = element_text(
               size = 12, family = mfont,
               color = "black"
          ),
          legend.margin = margin(t = -5, r = 0, b = 0, l = 0)
     )
print(f_starch)
ggsave("./Figs/4flour_vs_.png", f_starch,
     width = 5, height = 4,
     dpi = 130
)
#######
## Michaelis-Menten ###

# Estimation of kinetic parameters Vo, Km and Vmax
# ```{r kinetic, fig.align='center'}
# 1. Calculate Initial Velocity (v0)
# We calculate v0 using the change btw Time = 0 and Time = 15
v0_data <- UV2 %>%
     # Exclude the negative control
     filter(Starch != "Negative control") %>%
     # Filter only for the initial linear phase timepoints
     filter(Time %in% c(0, 15)) %>%
     # Group to average the two replicates in your dataset
     group_by(Starch, sub_conc, Time) %>%
     summarize(mean_ABS = mean(ABS), .groups = "drop") %>%
     # Reshape data to have Time 0 and Time 15 in separate columns
     pivot_wider(
          names_from = Time, values_from = mean_ABS,
          names_prefix = "T_"
     ) %>%
     mutate(
          # Calculate initial velocity
          v0 = (T_0 - T_15) / 15
     )

print(v0_data)

## Fitting Michaeli-Menten Model for the four types of starch
library(drc)
library(broom)
mm_results <- v0_data %>%
     group_by(Starch) %>%
     nest() %>%
     mutate(
          model = map(data, ~ nls(
               v0 ~ (Vmax * sub_conc) / (Km + sub_conc),
               data = .x,
               # Dynamic start values based on the data
               start = list(Vmax = max(.x$v0), Km = median(.x$sub_conc))
          )),
          # Extract the Vmax and Km estimates cleanly
          params = map(model, tidy)
     ) %>%
     unnest(params) %>%
     dplyr::select(Starch, term, estimate, std.error, p.value)

print(mm_results)

library(drc, quietly = TRUE)
fit_grouped <- drm(v0 ~ sub_conc,
     curveid = Starch,
     data = v0_data, fct = MM.2()
)
summary(fit_grouped)

pred_grid <- expand.grid(
     sub_conc = seq(0.1, 2.5, length.out = 100),
     Starch = unique(v0_data$Starch)
)
pred_grid$pred_v0 <- predict(fit_grouped, newdata = pred_grid)

######### plot the parameters side by side.
ggplot() +
     # Raw data points
     geom_point(data = v0_data, aes(
          x = sub_conc, y = v0,
          color = Starch
     ), size = 2.5) +
     # Fitted Michaelis-Menten lines
     geom_line(
          data = pred_grid, aes(
               x = sub_conc, y = pred_v0,
               color = Starch
          ),
          linewidth = 1
     ) +
     # Vertical dashed lines indicating Km
     geom_vline(
          data = mm_results |>
               filter(term == "Km"), aes(xintercept = estimate, color = Starch),
          linetype = "dashed", alpha = 0.7
     ) +
     geom_label(
          data = mm_results |>
               filter(term == "Km"), aes(x = estimate - c(0.2, 0, 0, 0), y = 0.01),
          label = "Km", size = 3, family = mfont, label.size = 0.2
     ) +
     facet_wrap(~Starch) +
     labs(
          x = "Substrate Concentration (%)",
          # y = expression("Initial Velocity " * v[0] * " (" * Delta * "Abs min"^"-1"*")")
          y = "Initial Velocity V<sub>0</sub> (&Delta;Abs min<sup>-1</sup>)"
     ) +
     theme_bw() +
     theme(
          text = element_text(
               family = "Times New Roman",
               color = "black",
               size = 12
          ),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 12, color = "black"),
          axis.title.y = element_markdown(color = "black"),
          legend.margin = margin(-10, 0, 0, 0),
          plot.margin = margin(5, 5, 5, 5),
          legend.position = "none",
          plot.caption = element_text(size = 12),
          strip.text = element_text(size = 12)
     )
ggsave("Figs/km_4starch.png",
     width = 5.5, height = 4,
     dpi = 300
)

########


facet_labels <- c(
     Km   = "K<sub>m</sub> (Substrate %)",
     Vmax = "V<sub>max</sub> (Abs min<sup>-1</sup>)"
)

mm_results |> ggplot(
     aes(Starch, estimate, fill = Starch)
) +
     geom_col(width = 0.6, alpha = 0.85, color = "black") +
     geom_errorbar(
          aes(ymin = estimate - std.error, ymax = estimate + std.error),
          width = 0.2
     ) +
     facet_wrap(. ~ term,
          scales = "free_y",
          labeller = labeller(
               term = as_labeller(facet_labels)
          )
     ) +
     labs(
          y = "",
          x = ""
     ) +
     theme_minimal() +
     theme(
          text = element_text(
               family = "Times New Roman",
               color = "black",
               size = 12
          ),
          # plot.margin = margin(5, 5, 5, 5),
          strip.text = element_markdown(
               family = "Times New Roman",
               size = 12,
               color = "black"
          ),
          axis.line = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 12),
          legend.position = "none",
          axis.text.x = element_text(angle = 50, hjust = 1),
          strip.background = element_rect(fill = "grey90")
     )

ggsave("Figs/km_vmax.png",
     width = 5.5, height = 4,
     dpi = 300
)
