library(tidyverse)
library(vegan)
library(FD)
library(picante)
library(rcompanion)
library(rstatix)
library(ggpubr)
library(ggpattern)


# load data
df <- read.csv("data/TA_Abundance.csv")

# Extracting environment variables
meta <- df[,1:4]
meta$region <- factor(meta$region, levels = c("Baikal","Penza", "Tobolsk", "Surgut", "Urengoi", "Zapolyarniy"))
meta$soil_type <- factor(meta$soil_type, levels = c("mineral soil","organic soil"))

# Extracting relative abundance table
cdm <- df[,5:181]
rowSums(cdm)
colnames(cdm) <- gsub("\\.", " ", colnames(cdm))
# cdm <- cdm[, sort(colnames(cdm))]


# Graphic format
plot.theme <- theme(plot.title = element_text(size = 16, color = "black", family = "sans", face= "bold"),
                    axis.line = element_line(size = .5, colour = "black"),
                    axis.ticks = element_line(color = "black"),
                    axis.text = element_text(size = 18, color = "black", family = "sans", vjust = 0.5, hjust = 0.5),
                    axis.title = element_text(size = 20, color = "black", family  = "sans", face = "bold", vjust = 0.5, hjust = 0.5),
                    legend.position = "none",
                    legend.background = element_blank(),
                    legend.key = element_blank(),
                    legend.text = element_text(colour = 'black', size = 20, family = "sans", face = 'bold'),
                    legend.title = element_blank(),
                    panel.background = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank())

#------------ Species diversity (SD) ------------#
S <- specnumber(cdm)
D1 <- exp(diversity(cdm))
D2 <- diversity(cdm, index = "invsimpson")
E <- diversity(cdm)/log(specnumber(cdm))

SD <- cbind.data.frame(S, D1, D2, E)


#------------Functional diversity and Community-weighted means of traits----------------#
# Load trait data
traits <- read.csv("data/TA_Traits.csv")
row.names(traits) <- traits$Species
row.names(traits) == colnames(cdm)
traits <- traits[,-1]

names(traits) <- c("L1","L2","L3","R1","R2","test.shape",
                   "ap1","ap2","R3","R4","ap.position","ap.inv",
                   "ap.rim","collar","cover","partition","spine","feed","Biovolume")

traits$test.shape <- factor(traits$test.shape, levels = c("sphere", "hemisphere", "cylinder", "patelliform", "rectangular cuboid", "ovoid", "pyriform", "spiral"))
traits$ap.position <- factor(traits$ap.position, levels = c("straight terminal", "sub terminal", "central ventral", "shifted ventral", "amphistomic"))
traits$ap.inv <- factor(traits$ap.inv, levels = c("absent", "slightly", "strongly"))
traits$ap.rim <- factor(traits$ap.rim, levels = c("straight", "curved", "lobbed", "denticular"))
traits$cover <- factor(traits$cover, levels = c("organic", "xenosomes", "idiosomes", "cleptosomes"))
traits$feed <- factor(traits$feed, levels = c("mixotrophy", "bacterivory", "predatory"))
traits$collar <- as.numeric(traits$collar == "presence")
traits$partition <- as.numeric(traits$partition == "presence")
traits$spine <- as.numeric(traits$spine == "presence")

# Calculate community-weighted means of biovolume
cwm <- as.matrix(cdm) %*% as.matrix(traits[, 19]) |> data.frame()
names(cwm) <- "CWM.biovolume" 

# Calculate functional diversity and Community-weighted means of traits
fd <- data.frame(dbFD(x = traits[, -19], a = as.matrix(cdm), w = c(rep(1/5, 5), 1, rep(1/8, 8), rep(1/3, 3), 1), CWM.type = "all"))

# Calculate SES.MPD
tr.d <- gowdis(traits[ ,-19], w = c(rep(1/5, 5), 1, rep(1/8, 8), rep(1/3, 3), 1))

ses.mpd <- ses.mpd(cdm, tr.d, null.model = "independentswap", runs = 999, iterations = 1000)
meta$mpd.z <- ses.mpd$mpd.obs.z

ses.mpd.a <- ses.mpd(cdm, tr.d, null.model = "independentswap", runs = 999, iterations = 1000, abundance.weighted = T)
meta$mpd.a.z <- ses.mpd.a$mpd.obs.z

# Aggregate all dependent variables
total <- cbind(meta[, c(1, 2)], SD, cwm, fd[, c(3,5,6,9:50)], meta[, c(5,6)])

median_result <- total |> 
  group_by(region, soil_type) |> 
  summarise(across(everything(), median, na.rm = TRUE)) |> 
  t() |> 
  data.frame()
# write.csv(median_result,"median_result.csv")

# Perform Scheirer-Ray-Hare test (non-parametric two-way ANOVA alternative)
# Tests main effects of soil type, region and their interaction on all dependent variables, get p-values
dv_names <- names(total)[3:54]

results <- data.frame(
  soil_type_p = character(length(dv_names)),
  region_p = character(length(dv_names)),
  interaction_p = character(length(dv_names)),
  row.names = dv_names,
  stringsAsFactors = FALSE
)

for (i in seq_along(dv_names)) {
  dv <- dv_names[i]
  current_formula <- formula(paste(dv, "~ soil_type + region"))
  
  tryCatch({

    srh_model <- scheirerRayHare(current_formula, data = total)
    p_vals_raw <- srh_model$p.value
    p_vals_formatted <- character(length(p_vals_raw))
    
    for (j in seq_along(p_vals_raw)) {
      if (is.na(p_vals_raw[j])) {
        p_vals_formatted[j] <- NA_character_
      } else if (p_vals_raw[j] < 0.001) {
        p_vals_formatted[j] <- "< 0.001"
      } else {
        p_vals_formatted[j] <- sprintf("%.3f", round(p_vals_raw[j], 3))
      }
    }
    
    if (length(p_vals_formatted) >= 3) {
      results[i, ] <- p_vals_formatted[1:3]
    } else {
      results[i, ] <- rep(NA_character_, 3)
    }
  }, error = function(e) {
    results[i, ] <<- rep(NA_character_, 3)
  })
}

print(results)
# write.csv(results,"p_results.csv")


#------------------------------visualization-----------------------------------#

# Define function to create boxplot with statistical annotations and data points
box.plot <- function(df, y) {
  # Perform Wilcoxon rank sum test between soil types within each region
  wt <- wilcox_test(group_by(df, region), as.formula(paste0(y, " ~ soil_type"))) 
  
  # Add significance labels to test results
  wt <- add_significance(wt, 'p')
  
  # Calculate positions for p-value annotations
  wt <- add_xy_position(wt, x = 'region')
  
  # Create base boxplot using ggpubr
  p <- ggboxplot(df, y = y, x = 'region', fill = 'soil_type', outlier.colour = NA) + 
    # Set custom fill colors for soil types
    scale_fill_manual(values = c('#FDB462','#6baed6')) +
    # Add jittered data points with semi-transparency
    geom_jitter(aes(fill = soil_type),   # Color points by soil type
                position = position_jitterdodge(  # Separate points by groups
                  jitter.width = 0.2,       # Horizontal jitter amount
                  dodge.width = 0.8),       # Align with boxplot dodging
                alpha = 0.6,                # 60% opacity
                size = 2,                   # Point size
                shape = 21,                 # Fillable circle shape
                color = "black") +          # Black border for points
    # Add statistical annotations
    stat_pvalue_manual(wt, 
                       color = "#d7301f",  # Red color for annotations
                       label = 'p.signif', # Show significance symbols
                       tip.length = 0.01,  # Length of annotation lines
                       size = 5) +         # Annotation symbol size
    plot.theme  
}


#-----------------------------Selected matrices for visualization

# Species richness (S)
# Perform Scheirer-Ray-Hare test (non-parametric two-way ANOVA alternative)
# Tests main effects of soil type, region and their interaction on species richness (S)
scheirerRayHare(S ~ soil_type + region, total)

# Record whether soil type, region and their interaction are significant or not, for graphing purposes
text <- paste0("Soil type: ", "p < 0.001", "; Region: ", "p < 0.001", "; Interaction: ", "p < 0.001")

# Fig. 3a
p <- box.plot(total, "S") + 
  theme(legend.position = c(0.85, 0.95)) + # Adjust legend position (x, y coordinates)
  labs(x = "Region", y = "Species richness", title = text) 
p

ggsave(p, filename = "figs/fig3a.png", width = 8, height = 6, dpi = 300)

# Other indices are annotated similar to species richness (S)

# Exponential of Shannon index D1
scheirerRayHare(D1 ~ soil_type + region, total)
text <- paste0("Soil type: ", "p < 0.001", "; Region: ", "p < 0.001", "; Interaction: ", "p < 0.001")

# Fig. 3b
p <- box.plot(total, "D1") + 
  labs(x = "Region", y = "Exponential of Shannon index", title = text)
p

ggsave(p, filename = "figs/fig3b.png", width = 8, height = 6, dpi = 300)

# Pielou evenness
scheirerRayHare(E ~ soil_type + region, total)
text <- paste0("Soil type: ", "p < 0.001", "; Region: ", "p < 0.001", "; Interaction: ", "p = 0.154")

# Fig. 3c
p <-  box.plot(total, "E") + 
  labs(x = "Region", y = "Pielou evenness", title = text)
p

ggsave(p, filename = "figs/fig3c.png", width = 8, height = 6, dpi = 300)

# Functional richness (FRic)
scheirerRayHare(FRic ~ soil_type + region, total)
text <- paste0("Soil type: ", "p < 0.001", "; Region: ", "p = 0.002", "; Interaction: ", "p < 0.001")

# Fig. 3d
p <- box.plot(total, "FRic") + 
  labs(x = "Region", title = text)
p

ggsave(p, filename = "figs/fig3d.png", width = 8, height = 6, dpi = 300)


#----------------------Community-weighted means (CWMs) for traits

# Biovolume
scheirerRayHare(CWM.biovolume ~ soil_type + region, total)
text <- paste0("Soil type: ", "p < 0.001", "; Region: ", "p < 0.001", "; Interaction: ", "p = 0.152")

# Fig. 4a
p <- box.plot(total, "CWM.biovolume") + 
  labs(x = "Region", y = "Biovolume (μm\u00b3)", title = text) + 
  theme(legend.position = c(0.7, 0.95)) + 
  ylim(0, 350000)
p

ggsave(p, filename = "figs/fig4a.png", width = 8, height = 6, dpi = 300)

# Shell width
scheirerRayHare(CWM.L2 ~ soil_type + region, total)
text <- paste0("Soil type: ", "p = 0.489", "; Region: ", "p < 0.001", "; Interaction: ", "p = 0.001")

# Fig. 4b
p <- box.plot(total, "CWM.L2") + 
  labs(x = "Region", y = "Shell width (μm)", title = text)
p

ggsave(p, filename = "figs/fig4b.png", width = 8, height = 6, dpi = 300)

# Ratio (shell width / shell length)
scheirerRayHare(CWM.R1 ~ soil_type + region, total) 
text <- paste0("Soil type: ", "p < 0.001", "; Region: ", "p = 0.001", "; Interaction: ", "p = 0.071")

# Fig. 4c
p <- box.plot(total, "CWM.R1") + 
  labs(x = "Region", y = "R1 (shell width / shell length)", title = text) 
p

ggsave(p, filename = "figs/fig4c.png", width = 8, height = 6, dpi = 300)

# Ovoid shell(or test) shape
scheirerRayHare(CWM.test.shape_ovoid ~ soil_type + region, total)
text <- paste0("Soil type: ", "p < 0.001", "; Region: ", "p < 0.001", "; Interaction: ", "p = 0.112")

# Fig. 4d
p <- box.plot(total, "CWM.test.shape_ovoid") + 
  labs(x = "Region", y = "Ovoid shell shape", title = text) 
p

ggsave(p, filename = "figs/fig4d.png", width = 8, height = 6, dpi = 300)


# Xenosomes shell covering
scheirerRayHare(CWM.cover_xenosomes ~ soil_type + region, total)
text <- paste0("Soil type: ", "p < 0.001", "; Region: ", "p < 0.001", "; Interaction: ", "p = 0.002")

# Fig. 4e
p <- box.plot(total,"CWM.cover_xenosomes") + 
  labs(x = "Region", y = "Xenosomes shell covering", title = text)
p

ggsave(p, filename = "figs/fig4e.png", width = 8, height = 6, dpi = 300)


# Aperture position as central ventral
scheirerRayHare(CWM.ap.position_central.ventral ~ soil_type + region, total)
text <- paste0("Soil type: ", "p < 0.001", "; Region: ", "p = 0.018", "; Interaction: ", "p = 0.056")

# Fig. 4f
p <- box.plot(total,"CWM.ap.position_central.ventral") + 
  labs(x = "Region", y = "Aperture position as central ventral", title =  text)
p

ggsave(p, filename = "figs/fig4f.png", width = 8, height = 6, dpi = 300)


# Denticular aperture rim
scheirerRayHare(CWM.ap.rim_denticular ~ soil_type + region, total)
text <- paste0("Soil type: ", "p = 0.003", "; Region: ", "p < 0.001", "; Interaction: ", "p = 0.634")

# Fig. 4g
p <- box.plot(total,"CWM.ap.rim_denticular") + 
  labs(x = "Region", y = "Denticular aperture rim", title = text)
p

ggsave(p, filename = "figs/fig4g.png", width = 8, height = 6, dpi = 300)



# Feeding type as predator
scheirerRayHare(CWM.feed_predatory ~ soil_type + region, total)
text <- paste0("Soil type: ", "p < 0.001", "; Region: ", "p = 0.001", "; Interaction: ", "p < 0.001")

# Fig. 4h
p <- box.plot(total,"CWM.feed_predatory") + 
  labs(x = "Region", y = "Feeding type as predator", title = text)
p

ggsave(p, filename = "figs/fig4h.png", width = 8, height = 6, dpi = 300)


# SES.MPD
scheirerRayHare(mpd.z ~ soil_type + region, meta) # Non-weighted SES.MPD
scheirerRayHare(mpd.a.z ~ soil_type + region, meta) # Weighted SES.MPD

meta.med <- meta |> group_by(soil_type, region) |> 
  summarise(ses.mpd = median(mpd.z, na.rm = T), 
            pp = ifelse(wilcox.test(mpd.z)$p.value < 0.05, "non-random", "random"),
            ses.mpd.a = median(mpd.a.z, na.rm = T),
            pp.a = ifelse(wilcox.test(mpd.a.z)$p.value < 0.05, "non-random", "random"))

meta$id <- paste0(meta$soil_type, "-", meta$region)
meta.med$id <- paste0(meta.med$soil_type, "-", meta.med$region)
meta.join <- meta |> select(id) |> left_join(meta.med)
df.null <- cbind(meta.join[, c(2:3,5,7)], meta[, 5:6])

# Non-weighted SES.MPD
# Fig. 5a
text <- paste0("Soil type: ", "p < 0.001", "; Region: ", "p < 0.001", "; Interaction: ", "p < 0.001")
p <- ggplot(df.null, aes(x = region, y = mpd.z, 
                         pattern = pp, fill = soil_type, 
                         group = interaction(region, soil_type))) +  
  geom_boxplot(
    aes(fill = soil_type),
    position = position_dodge(preserve = "single", width = 0.8),
    width = 0.7,
    outlier.shape = NA  
  ) +
  geom_boxplot_pattern(
    aes(pattern = pp, pattern_fill = pp),
    position = position_dodge(preserve = "single", width = 0.8),
    color = "black", 
    pattern_fill = "white", 
    pattern_angle = 45, 
    pattern_density = 0.1, 
    pattern_spacing = 0.025, 
    pattern_key_scale_factor = 0.6,
    width = 0.7,
    outlier.shape = NA  
  ) +
  geom_jitter(
    aes(fill = soil_type),  
    position = position_jitterdodge(  
      dodge.width = 0.8,     
      jitter.width = 0.2,    
      jitter.height = 0      
    ),
    shape = 21,             
    color = "black",       
    size = 2,            
    alpha = 0.6           
  ) +
  scale_fill_manual(name = "soil_type", values = c('#FDB462', '#80B1D3')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) + 
  guides(
    pattern = guide_legend(override.aes = list(fill = "white")),
    fill = guide_legend(override.aes = list(pattern = "none"))
  ) +
  plot.theme +
  theme(legend.position=c(0.88,0.88)) +
  labs(x = "Region", y = "SES.MPD", title = text)
p

ggsave(p, filename = "figs/fig5a.png", width = 8, height = 6, dpi = 300)


# Weighted SES.MPD
# Fig. 5b
text <- paste0("Soil type: ", "p < 0.001", "; Region: ", "p < 0.001", "; Interaction: ", "p < 0.001")
p <- ggplot(df.null, aes(x = region, y = mpd.a.z, 
                         pattern = pp.a, fill = soil_type, 
                         group = interaction(region, soil_type))) +  
  geom_boxplot(
    aes(fill = soil_type),
    position = position_dodge(preserve = "single", width = 0.8),
    width = 0.7,
    outlier.shape = NA  
  ) +
  geom_boxplot_pattern(
    aes(pattern = pp.a, pattern_fill = pp.a),
    position = position_dodge(preserve = "single", width = 0.8),
    color = "black", 
    pattern_fill = "white", 
    pattern_angle = 45, 
    pattern_density = 0.1, 
    pattern_spacing = 0.025, 
    pattern_key_scale_factor = 0.6,
    width = 0.7,
    outlier.shape = NA  
  ) +
  geom_jitter(
    aes(fill = soil_type),  
    position = position_jitterdodge(  
      dodge.width = 0.8,     
      jitter.width = 0.2,    
      jitter.height = 0      
    ),
    shape = 21,             
    color = "black",       
    size = 2,            
    alpha = 0.6           
  ) +
  scale_fill_manual(name = "soil_type", values = c('#FDB462', '#80B1D3')) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) + 
  guides(
    pattern = guide_legend(override.aes = list(fill = "white")),
    fill = guide_legend(override.aes = list(pattern = "none"))
  ) +
  plot.theme +
  labs(x = "Region", y = "SES.MPD (abundance-weighted)", title = text)
p

ggsave(p, filename = "figs/fig5b.png", width = 8, height = 6, dpi = 300)
