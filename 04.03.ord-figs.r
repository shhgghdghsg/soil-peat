library(tidyverse)
library(vegan)
library(pairwiseAdonis)

# load data
df <- read.csv("data/TA_Abundance.csv")

# Extracting environment variables
meta <- df[, 1:4]
meta$region <- factor(meta$region, levels = c("Baikal","Penza", "Tobolsk", "Surgut", "Urengoi", "Zapolyarniy"))
meta$soil_type <- factor(meta$soil_type, levels = c("mineral soil", "organic soil"))

meta$id <- paste0(meta$soil_type, sep = "-", meta$region)
meta$id <- factor(meta$id, levels = c("organic soil-Baikal","organic soil-Penza","organic soil-Tobolsk","organic soil-Surgut","organic soil-Urengoi","organic soil-Zapolyarniy",
                                      "mineral soil-Baikal","mineral soil-Penza","mineral soil-Tobolsk","mineral soil-Surgut","mineral soil-Urengoi","mineral soil-Zapolyarniy"))

# Extracting relative abundance table
cdm <- df[,5:181]
rowSums(cdm)
colnames(cdm) <- gsub("\\.", " ", colnames(cdm))


# PERMANOVA
panova <- adonis2(cdm ~ soil_type * region, meta, method = "bray", permutations = 9999)
panova

# Pairwise comparison
pair.res <- pairwise.adonis2(cdm ~ id, meta, method = "bray", permutations = 9999)
pair.res
# write.csv(pair.res, "pairwise.csv")


# Species Hellinger Transformation
spe.h <- decostand(cdm, "hellinger")

# PCA
spe.h.pca <- rda(spe.h)

(axis1 <- round(spe.h.pca$CA$eig[1]/sum(spe.h.pca$CA$eig)*100,2))
(axis2 <- round(spe.h.pca$CA$eig[2]/sum(spe.h.pca$CA$eig)*100,2))
axis1 + axis2

sp.sc <- scores(spe.h.pca, display = "species") 
st.sc <- scores(spe.h.pca, display = "sites")

st.sc <- cbind(as.data.frame(st.sc), zone = meta$id)

cent <- aggregate(cbind(PC1, PC2) ~ zone, data = st.sc, FUN = mean)
cent$text <- c("OB","OP","OT","OS","OU","OZ", "MB","MP","MT","MS","MU","MZ")
cent$pos1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, -0.005, 0, 0)
cent$pos2 <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.015, 0)

segs <- merge(st.sc, setNames(cent, c('zone','oPC1','oPC2')),
              by = 'zone', sort = FALSE)


# Create a dataframe to extract the convex hull points
pca_hull <- 
  st.sc |> 
  group_by(zone) |> 
  slice(chull(PC1, PC2))


# Graphic format
plot.theme <- theme(plot.title = element_text(size = 16, color = "black", family = "sans", face = "bold"),
                   axis.line = element_line(linewidth = .5, colour = "black"),
                   axis.ticks = element_line(color = "black"),
                   axis.text = element_text(size = 18, color = "black", family = "sans", vjust = 0.5, hjust = 0.5),
                   axis.title = element_text(size = 20, color = "black", family = "sans",face = "bold", vjust = 0.5, hjust = 0.5),
                   legend.position = c(0.15, 0.81),
                   legend.background = element_blank(),
                   legend.key = element_blank(),
                   legend.text = element_text(colour = 'black', size = 20, family = "sans", face = 'bold'),
                   legend.title = element_blank(),
                   panel.background = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())


text <- bquote('PERMANOVA (Soil type: p < 0.001,' ~~ R^'2' ~~ '= 9.4%; ' ~~ 'Region: p < 0.001, ' ~~ R^'2' ~~ '= 13.6%; ' ~~ 'Interaction: p < 0.001,' ~~ R^'2' ~~ '= 10.0%)')

# Fig. 2a
p <- ggplot(st.sc, aes(x = PC1, y = PC2, colour = zone)) +
  geom_segment(data = segs, mapping = aes(xend = oPC1, yend = oPC2)) + # spiders
  geom_point(data = cent, size = 5) +
  geom_point() +
  geom_hline(yintercept = 0, lty = 2,col = "black") + 
  geom_vline(xintercept = 0, lty = 2,col = "black") +
  labs(x = "PC1(16.5%)", y = "PC2(10.9%)", title = text) +
  annotate('text', x = cent$PC1 + cent$pos1, y = cent$PC2 + cent$pos2, label = cent$text, size = 6) +
  plot.theme +
  guides(colour = guide_legend(override.aes = list(shape = 15)))+
  scale_color_discrete(labels = c("OB: organic soil-Baikal",
                                  "OP: organic soil-Penza",
                                  "OT: organic soil-Tobolsk",
                                  "OS: organic soil-Surgut",
                                  "OU: organic soil-Urengoi",
                                  "OZ: organic soil-Zapolyarniy",
                                  "MB: mineral soil-Baikal",
                                  "MP: mineral soil-Penza",
                                  "MT: mineral soil-Tobolsk",
                                  "MS: mineral soil-Surgut",
                                  "MU: mineral soil-Urengoi",
                                  "MZ: mineral soil-Zapolyarniy"))

p
ggsave(p, filename = "figs/fig2a.png", width = 15, height = 10, dpi = 600)


# Fifteen species with 15 most influential species (with the longest arrows) 
len <- sqrt(rowSums(sp.sc^2))
top <- order(len, decreasing = T)[1:15]

sp.sc.top <- data.frame(sp.sc[top, ])

sp.sc.top$vjust <- c(0.15,-0.02,0,0,0,0.14,0.04,-0.05,-0.12,-0.12,0.13,0.095,0.12,0,-0.13)
sp.sc.top$hjust <- c(0.03,-0.03,-0.02,0.03,-0.02,0,-0.03,-0.03,-0.01,0.03,0,-0.01,0.01,-0.02,0)

# Fig. 2b
p <- ggplot(st.sc, aes(x = PC1, y = PC2, colour = zone)) +
  geom_segment(data = sp.sc.top, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(angle = 22.5, length = unit(0.2, "cm"), type = "closed"), 
               linetype = 1, linewidth = 0.6, color = "black") +
  geom_hline(yintercept = 0, lty = 2,col = "black") + 
  geom_vline(xintercept = 0, lty = 2,col = "black") +
  labs(x = "PC1(16.5%)", y = "PC2(10.9%)") +
  annotate('text', x = sp.sc.top$PC1 + sp.sc.top$vjust, y = sp.sc.top$PC2 + sp.sc.top$hjust, label = rownames(sp.sc.top), size = 6, color = "blue") +
  plot.theme +
  theme(legend.position="none")

p
ggsave(p, filename = "figs/fig2b.png", width = 15, height = 10, dpi = 600)