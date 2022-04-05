

library(phyloseq)
library(ampvis2)
library(vegan)
library(microbiome)
library(microbiomeutilities)
library(RColorBrewer)
library(metagMisc)
library("ape")
library(dplyr)
library(ggpubr)
library(ggstatsplot)
library(ggpubr)
library(patchwork)
library(lsmeans)
library(agricolae)
library(rstatix)
library(multcompView)
library(stringr)

devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")

# centrifuge

phylo = readRDS("./phylo_original.rds")

phylo #96 samples

meta = read.csv("./metadata_soil.csv", row.names = 1)
sample_data(phylo) = meta

# # check data
min(taxa_sums(phylo)) #1
min(sample_sums(phylo)) # 408
# what samples would be lost
sum(sample_sums(phylo) >= 20000) #92 at 20000
sample_names(phylo)[sample_sums(phylo) < 3000]

# plot rarefaction curve
amp.phylo = phyloseq_to_ampvis2(phylo)
rarecu = ampvis2::amp_rarecurve(amp.phylo)
rarecu
rarecu + geom_vline(aes(xintercept=20000), color="red")
rarecu + geom_vline(aes(xintercept=30000), color="green")


phylo.subset = subset_taxa(phylo, !is.na(Kingdom) & !is.na(Phylum))
phylo.trim = subset_taxa(phylo.subset, Phylum != "Cyanobacteria")

phylo.prune0 = prune_taxa(taxa_sums(phylo.trim) > 10, phylo.trim)
phylo.prune0 #9740
phylo.prune = prune_samples(sample_sums(phylo.prune0) >= 20000, phylo.prune0)
phylo.prune

# save prune data
#saveRDS(phylo.prune, "./phylo_prune.rds")

# remove Cyanobacteria
tax_table(phylo.prune)[, "Kingdom"]

# filter out other and keep bacteria only
phylo.bac = subset_taxa(phylo.prune, Kingdom == "Bacteria")
phylo.bac #9312

# saveRDS(phylo.bac, "./phylo_prune_tax10_sam20000_bac.rds")

phylo.bac = readRDS("./phylo_prune_tax10_sam20000_bac.rds")





#-----------------------------------------------------------------------------------------

psq2 = subset_samples(phylo.bac, Year=="2017" | Year=="2018")

# format to best hit
bst_psq2 = format_to_besthit(psq2)

# saveRDS(psq2, "./psq2_bac.rds")
# saveRDS(bst_psq2, "./bst_psq2_bac.rds")

psq2 = readRDS("./psq2_bac.rds")


psq2
psq2.melt = psq2 %>%
  psmelt()

head(psq2.melt)




psq2.melt %>%
  group_by(Year, Season) %>%
  summarise(ndist = n_distinct(Genus))

summarize_phyloseq(psq2)
mean(readcount(psq2))

top_taxa(psq2, n=10)

amp.psq2 = phyloseq_to_ampvis2(psq2)

amp.psq2

ampRARE = amp_rarecurve(amp.psq2, color_by = "State", facet_by = "Season")
ampRARE +
  theme(text = element_text(size=15),
        axis.title.x = element_text(margin = margin(t=20)),
        axis.title.y = element_text(margin = margin(r=20)))



sum(sample_sums(psq2)) #7519765
meta2 = data.frame(sample_data(psq2))

head(meta2)
min(sample_sums(psq2))

min(sample_sums(psq2))

# rarefaction
rare = rarefy_even_depth(psq2, rngseed = 123, sample.size = min(sample_sums(psq2)), replace = F)
rch = estimate_richness(rare, measures = c("Observed", "Shannon", "InvSimpson"))

meta_al = merge(meta2, rch, by = "row.names")
head(meta_al)



# gather the columns
adiv = meta_al %>%
  tidyr::gather(key="indices", value="diversity", Observed:InvSimpson)

# saveRDS(adiv, file = "./gathered_alpha_diversity.rds")

# rm(list = ls())

rm(list = ls())
adiv = readRDS("./gathered_alpha_diversity.rds")
head(adiv)


#--------------------------------------------------------------------------------
# alpha diversity statistics
#--------------------------------------------------------------------------------
head(adiv)
adiv %>%
  dplyr::filter(indices == "Observed") %>%
  dplyr::group_by(Year, State, Season) %>%
  dplyr::summarise(avg = mean(diversity), std = sd(diversity)/sqrt(n())) %>% 
  write.csv("./output/alpha_diversity/Observed_groups_bac.csv", quote = F, row.names = F)

adiv %>%
  dplyr::filter(indices == "Shannon") %>%
  dplyr::group_by(Year, State, Season) %>%
  dplyr::summarise(avg = mean(diversity), std = sd(diversity)/sqrt(n())) %>% 
  write.csv("./output/alpha_diversity/Shannon_groups_bac.csv", quote = F, row.names = F)


source("./soilmicrobiome_helpers.R")


adiv$Year = as.factor(adiv$Year)
adiv$Season = as.factor(adiv$Season)
adiv$State = as.factor(adiv$State)

# observed
obs_adiv2017 = adiv %>%
  filter(indices == "Observed" & Year == "2017")

#
obs_adiv_out2017.2 = aov(diversity ~ State + Season + Season : State, data = obs_adiv2017)
summary(obs_adiv_out2017.2)

# tukey test and plot
(obs2017.2 = plotTukeyDive(obs_adiv2017, obs_adiv_out2017.2, ylabs = "Observed OTU", facet_by = "~ Season", LetterTerm = "State:Season", plotx = "State", groupbys = c("State", "Season"), xlabs = ""))

getDiff(obs_adiv2017)

#
obs_adiv2018 = adiv %>%
  filter(indices == "Observed" & Year == "2018")


obs_adiv_out2018.2 = aov(diversity ~ State + Season + Season : State, data = obs_adiv2018)
summary(obs_adiv_out2018.2)

#
(obs2018.2 = plotTukeyDive(obs_adiv2018, obs_adiv_out2018.2, ylabs = "Observed OTU", facet_by = "~ Season", LetterTerm = "State:Season", plotx = "State", groupbys = c("State", "Season")))

ggarrange(obs2017.2$fig, obs2018.2$fig, labels = c("A", "B"), ncol = 1, align = "hv")

getDiff(obs_adiv2018)

# shannon index
(shan_adiv2017 = adiv %>%
    filter(indices == "Shannon" & Year == "2017"))

shan_aov_out2017.2 = aov(diversity ~ State + Season + Season : State, data = shan_adiv2017)
summary(shan_aov_out2017.2)

(shan2017.2 = plotTukeyDive(shan_adiv2017, shan_aov_out2017.2, gap=0.1, ylabs = "Shannon index", facet_by = "~ Season", LetterTerm = "State:Season", plotx = "State", groupbys = c("State", "Season")))

getDiff(shan_adiv2017)

(shan_adiv2018 = adiv %>%
    filter(indices == "Shannon" & Year == "2018"))


#
shan_aov_out2018.2 = aov(diversity ~ State + Season + Season : State, data = shan_adiv2018)
summary(shan_aov_out2018.2)

(shan2018.2 = plotTukeyDive(shan_adiv2018, shan_aov_out2018.2, gap = 0.1, ylabs = "Shannon index", facet_by = "~ Season", LetterTerm = "State:Season", plotx = "State", groupbys = c("State", "Season")))

ggarrange(shan2017.2$fig, shan2018.2$fig, labels = c("A", "B"), ncol = 1)

getDiff(shan_adiv2018)


#-----------------------------------------------------------------------------------
# beta diversity
#-----------------------------------------------------------------------------------


# beta diversity
rm(list = ls())
bst_psq2 = readRDS("./bst_psq2_bac.rds")

# subset samples to 2017 and 2018
psq2017 = subset_samples(bst_psq2, Year=="2017")
psq2018 = subset_samples(bst_psq2, Year=="2018")

head(sample_data(psq2017))

DistBC2017 = phyloseq::distance(psq2017, method = "bray")
ordBC2017 = ordinate(psq2017, method = "PCoA", distance = DistBC2017)

DistBC2018 = phyloseq::distance(psq2018, method = "bray")
ordBC2018 = ordinate(psq2018, method = "PCoA", distance = DistBC2018)


jsdf2017 = plot_ordination(psq2017, ordBC2017, color = "State", shape = "Season",justDF = T)
jsdf2018 = plot_ordination(psq2018, ordBC2018, color = "State", shape = "Season", justDF = T)

# ploting
(beta2017 = jsdf2017 %>% ggplot(aes(x=Axis.1, y=Axis.2, color=Season, shape=State, group = Season)) +
    geom_point(size=6, alpha=0.5) +
    theme_classic() +
    stat_ellipse(alpha=0.5) +
    # geom_polygon(data=hull2017, color="grey", alpha=0.02) +
    # facet_grid(~Year) +
    labs(color="Season", shape="State", x=paste0("Axis.1 [ ", round(100*ordBC2017$values[1, 2],1) ,"% ]"), y=paste0("Axis.2 [ ", round(100*ordBC2017$values[2, 2],1) ,"% ]"))+
    theme(text=element_text(size=12),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.title.x = element_text(margin = margin(t=5, r=0, b=0, l=0)),
          axis.title.y = element_text(margin = margin(t=0, r=5, b=0, l=0))))

(beta2018 = jsdf2018 %>% ggplot(aes(x=Axis.1, y=Axis.2, color=Season, shape=State, group = Season)) +
    geom_point(size=6, alpha=0.5) +
    theme_classic() +
    stat_ellipse(alpha=0.5) +
    # geom_polygon(data=hull2018, alpha=0.05, color="grey") +
    # facet_grid(~Year) +
    labs(color="Season", shape="State", x=paste0("Axis.1 [ ", round(100*ordBC2018$values[1, 2],1) ,"% ]"), y=paste0("Axis.2 [ ", round(100*ordBC2018$values[2, 2],1) ,"% ]")) +
    theme(text=element_text(size=12),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          axis.title.x = element_text(margin = margin(t=5, r=0, b=0, l=0)),
          axis.title.y = element_text(margin = margin(t=0, r=5, b=0, l=0))))

ggarrange(beta2017, beta2018, labels=c("A", "B"),  common.legend = T, legend = "right")

# statistical analysis of betadiversity using adonis and disper
meta_al = data.frame(sample_data(bst_psq2))

meta_al2017 = meta_al %>%
  filter(Year=="2017")
head(meta_al2017)

meta_al2018 = meta_al %>%
  filter(Year=="2018")
head(meta_al2018)


# dispersion 2017
permdist2017.season= betadisper(DistBC2017, group = meta_al2017$Season)
anova(permdist2017.season)
plot(permdist2017.season)
boxplot(permdist2017.season$distances~permdist2017.season$group)

permdist2017.state= betadisper(DistBC2017, group = meta_al2017$State)
anova(permdist2017.state)
plot(permdist2017.state)
boxplot(permdist2017.state$distances~permdist2017.state$group)

# 2018
permdist2018.season= betadisper(DistBC2018, group = meta_al2018$Season)
anova(permdist2018.season)
plot(permdist2018.season)
boxplot(permdist2018.season$distances~permdist2018.season$group)

permdist2018.state= betadisper(DistBC2018, group = meta_al2018$State)
anova(permdist2018.state)
plot(permdist2018.state)
boxplot(permdist2018.state$distances~permdist2018.state$group)


perms = 10000
set.seed(10001)
perm2017.season = adonis(DistBC2017~Season+State, data=meta_al2017, permutations = perms)
perm2017.season
perm2017.season.2 = adonis(DistBC2017~State+Season, data=meta_al2017, permutations = perms)
perm2017.season.2


perm2018.season = adonis(DistBC2018~Season+State, data=meta_al2018, permutations = perms)
perm2018.season
perm2018.season.2 = adonis(DistBC2018~State+Season, data=meta_al2018, permutations = perms)
perm2018.season.2

###############################################################################
# dendrogram
# beta diversity
rm(list = ls())

bst_psq2 = readRDS("./bst_psq2_bac.rds")

psq2017.es = subset_samples(bst_psq2, Year== "2017" & Season == "ES")
psq2017.lf = subset_samples(bst_psq2, Year=="2017" & Season == "LF")

psq2018.es = subset_samples(bst_psq2, Year=="2018" & Season == "ES")
psq2018.lf = subset_samples(bst_psq2, Year=="2018" & Season == "LF")


source("./soilmicrobiome_helpers.R")
d2017.es = plotDendro(psq2017.es)
d2017.lf = plotDendro(psq2017.lf, yLab = "")
d2018.es = plotDendro(psq2018.es, xLab = "State")
d2018.lf = plotDendro(psq2018.lf, yLab = "", xLab = "State")

ggarrange(d2017.es, d2017.lf, d2018.es, d2018.lf, labels = c("A", "B", "C", "D"), align = "hv")


######################################################
source("./soilmicrobiome_helpers.R")


(bray2017.season = plotBray(psq2017, "Season", textfontsize = 5, gap = 0.05, xlabs = "Season"))
(bray2017.state = plotBray(psq2017, "State", textfontsize = 5, gap = 0.05, xlabs = "State"))
(bray2018.season = plotBray(psq2018, "Season", textfontsize = 5, ylabs="", gap = 0.05, xlabs="Season"))
(bray2018.state = plotBray(psq2018, "State", textfontsize = 5, ylabs = "", gap = 0.05, xlabs = "State"))

ggarrange(bray2017.season, bray2018.season, bray2017.state, bray2018.state, labels = c("A", "B", "C", "D"), legend = "none")

#-----------------------------------------------------------------------------------

# composition

#-----------------------------------------------------------------------------------
rm(list = ls())

bst_psq2 = readRDS("./bst_psq2_bac.rds")
taxlevel = "Order"

# unique_orders = tax_table(bst_psq2)[,"Order"]
# unique_orders.df = as.data.frame(unique_orders)
# uni_ord = unique(unique_orders.df$Order)


melt = bst_psq2 %>%
  tax_glom(taxlevel) %>%
  psmelt()

#melt %>%
#  write.csv("./output/ver3/relative_abundance_order_level_all.csv", quote = F, row.names = F)

perc = 1
thres = paste0("< ", perc, "%")

(overal_top10 = melt %>%
    dplyr::mutate(total_abund = sum(Abundance)) %>%
    dplyr::select(Season, State, !!sym(taxlevel), Abundance, total_abund) %>%
    dplyr::group_by(!!sym(taxlevel)) %>%
    dplyr::mutate(rel_abund = 100*Abundance / total_abund) %>%
    dplyr::mutate(genus_lab = ifelse(rel_abund < perc, thres, as.character(!!sym(taxlevel)))) %>%
    dplyr::summarise(sum_abund = sum(rel_abund)) %>%
    dplyr::arrange(desc(sum_abund)))

overal_top10 %>%
  write.table("./output/ver4/overall_top10_relative_abundance.tsv", sep = "\t", quote = F, row.names = F)

source("./soilmicrobiome_helpers.R")

overal_top10$Order = as.character(overal_top10$Order)
str(overal_top10)
head(overal_top10, n=10)

overal.top10.2017= kw_test(melt, bac_vector = overal_top10$Order[1:10], out = T, outdir = "./output/ver4/")
overal.top10.2018= kw_test(melt, bac_vector = overal_top10$Order[1:10], out = T, outdir = "./output/ver4/", filterYear = 2018)

overal.top10.2017$testR %>%
  write.csv("./output/ver4/top10_orders_testKW.csv", quote = F, row.names = F)

#order_rd %>%
#  write.table("./output/ver3/relative_abundance_at_order_level.tsv", quote = F, row.names = F, sep = "\t")


# melt %>%
#   dplyr::mutate(total_abund = sum(Abundance), rel_abund = 100*Abundance / total_abund) %>%
#   dplyr::select(Season, State, !!sym(taxlevel), Abundance, total_abund, rel_abund) %>%
#   dplyr::mutate(genus_lab = ifelse(rel_abund < perc, thres, as.character(!!sym(taxlevel)))) %>%
#   group_by(Order) %>%
#   dplyr::summarise(sum_abund = sum(rel_abund)) %>%
#   dplyr::arrange(desc(sum_abund)) %>%
#   write.csv("./output/ver3/relative_abundance_at_order_level.csv", quote = F)
#   


(top10.ord.df = melt %>%
    dplyr::filter( Domain == "Bacteria") %>%
    dplyr::group_by(Year, Season, State) %>%
    dplyr::mutate(total_abund = sum(Abundance)) %>%
    dplyr::select(Season, State, !!sym(taxlevel), Abundance, total_abund) %>%
    dplyr::group_by(Year, Season, State, !!sym(taxlevel)) %>%
    dplyr::mutate(rel_abund = 100*Abundance / total_abund) %>%
    dplyr::mutate(genus_lab = ifelse(rel_abund < perc, thres, as.character(!!sym(taxlevel))))  %>%
    dplyr::summarise(relative_abundance = sum(rel_abund), mean_relative_abundance = mean(rel_abund))%>%
    slice_max(order_by = relative_abundance, n=10))

top10.ord.df %>%
  write.table("./output/ver4/top10_orders.tsv",sep="\t", quote = F, row.names = F)

(top10.ord.df.byyear = melt %>%
    dplyr::filter( Domain == "Bacteria") %>%
    dplyr::group_by(Year, Season) %>%
    dplyr::mutate(total_abund = sum(Abundance)) %>%
    dplyr::select(Season, State, !!sym(taxlevel), Abundance, total_abund) %>%
    # mutate(genus_lab = ifelse(rel_abund < perc, thres, as.character(!!sym(taxlevel))))  %>%
    dplyr::group_by(Year, Season, !!sym(taxlevel)) %>%
    mutate(rel_abund = 100*Abundance / total_abund) %>%
    dplyr::summarise(relative_abundance = sum(rel_abund), mean_relative_abundance = mean(rel_abund))%>%
    slice_max(order_by = relative_abundance, n=10) %>%
    write.table("./output/ver4/top_order_byYearSeason.tsv", sep="\t",quote = F, row.names = F))
top10.ord.df.byyear


# this is to make the plot
(top6.ord.df = melt %>%
    dplyr::filter( Domain == "Bacteria") %>%
    dplyr::group_by(Year, Season, State) %>%
    dplyr::mutate(total_abund = sum(Abundance)) %>%
    dplyr::select(Season, State, !!sym(taxlevel), Abundance, total_abund) %>%
    dplyr::group_by(Year, Season, State, !!sym(taxlevel)) %>%
    dplyr::mutate(rel_abund = 100*Abundance / total_abund) %>%
    mutate(genus_lab = ifelse(rel_abund < perc, thres, as.character(!!sym(taxlevel)))) %>%  
    dplyr::group_by(Year, Season, State, genus_lab) %>%
    dplyr::summarise(relative_abundance = sum(rel_abund))%>%
    filter(genus_lab != thres) %>%
    slice_max(order_by = relative_abundance, n=3))


top6.ord = unique(top6.ord.df$genus_lab)[1:6]
top6.ord

source("./soilmicrobiome_helpers.R")


top6.ord.2 = factor(top6.ord)


# melt %>%
#   dplyr::group_by(Year, Season, State) %>%
#   dplyr::mutate(total_abund = sum(Abundance), rel_abund = 100*Abundance / total_abund) %>%
#   dplyr::select(Season, State, !!sym(taxlevel), Abundance, total_abund, rel_abund) %>%
#   mutate(genus_lab = ifelse(rel_abund < perc, thres, as.character(!!sym(taxlevel)))) %>%
#   dplyr::group_by(Year, Season, State, genus_lab) %>%
#   dplyr::summarise(relative_abundance = sum(rel_abund)) %>%
#   filter(Year==2017 & State == "NY" & genus_lab == "Burkholderiales")


kw2017 = kw_test(melt, bac_vector = top6.ord.2, out = F, outdir = "./output/ver4/kw_test/")

# df = kw2017$data
# 
# df$Burkholderiales %>%
#   filter(State == "NY") %>%
#   dplyr::group_by(Season) %>%
#   dplyr::summarise(mean_abund = sum(rel_abund))

figurelist2017 = plot_kw(kw2017, top6.ord.2, ylab = "", annoX = 0.5,xlab = "")

fig.top.2017 = ggarrange(plotlist = figurelist2017$figlist, align = "hv", common.legend = T, legend = "right") 
kw2017_A = annotate_figure(fig.top.2017, bottom = text_grob("State"),left = text_grob("Relative abundance (%)",  rot = 90, size=12))
kw2017_A

kw2018 = kw_test(melt, bac_vector = top6.ord.2, filterYear = 2018,out = F, outdir = "./output/ver4/kw_test/")

figurelist2018 = plot_kw(kw2018, top6.ord.2, ylab = "", annoX = 0.5, xlab = "")
fig.top.2018 = ggarrange(plotlist = figurelist2018$figlist, align = "hv", common.legend = T, legend = "right")
kw2018_B = annotate_figure(fig.top.2018, bottom = text_grob("State"), left = text_grob("Relative abundance (%)",  rot = 90, size=12))
kw2018_B

ggarrange(kw2017_A, kw2018_B, labels = c("A", "B"), align = "hv",nrow = 2, common.legend = T, legend = "right")
# ggarrange(kw2018_B, labels = "B")

#-----------------------------------------------------------------
source("./soilmicrobiome_helpers.R")


top10.ord = unique(top10.ord.df$genus_lab)
top10.ord.2 = factor(top10.ord)
top10.ord.2 = top10.ord.2[top10.ord.2 != "< 1%"]
top10.ord.2

kw2017.top = kw_test(melt, bac_vector = top10.ord.2, out = F, outdir = "./output/ver4/kw_test/")

kw2017.top

df2017.top = aggreDF(melt, yearfilter = 2017, taxfilter = "Bacteria", ntop = 10, groups=c("Season", "State"), perc = 1, taxlevel = taxlevel, sigStart = T, kw_out = kw2017.top)

unique(df2017.top$A$genus_lab)

kw2018.top = kw_test(melt, bac_vector = top10.ord.2, filterYear = 2018, out = F, outdir = "./output/kw_test/")

df2018.top = aggreDF(melt, yearfilter = 2018, taxfilter = "Bacteria", ntop = 10, groups=c("Season", "State"), perc = 1, taxlevel = taxlevel, sigStart = T, kw_out = kw2018.top)

unique(df2018.top$A$genus_lab)

df2018.top$A[df2018.top$A$Order == "Bacillales",]


colorIN2017 = colorPats(df2017.top,special = df2017.top$B)
colorIN2017

colorIN2018 = colorPats(df2018.top,special = df2018.top$B)
colorIN2018

# add * to top order in the legend, */* both season and state, */ season /* state

(es2017 = df2017.top$A %>%
    filter(Season == "ES") %>%
    plotComp(., plotx = "State", taxlevel = taxlevel, colorPalette = colorIN2017, facet = "", fill_legend = taxlevel, plot_title = "Early summer"))

(lf2017 = df2017.top$A %>%
    filter(Season == "LF") %>%
    plotComp(., plotx = "State", taxlevel = taxlevel, colorPalette = colorIN2017, facet = "", fill_legend = taxlevel, plot_title = "Late fall", ylabs = ""))

(es2018 = df2018.top$A %>%
    filter(Season == "ES") %>%
    plotComp(., taxlevel = taxlevel, colorPalette = colorIN2018, plotx = "State", facet = "",
             fill_legend = taxlevel, plot_title = "", xLabs = "State"))

(lf2018 = df2018.top$A %>%
    filter(Season == "LF") %>%
    plotComp(., taxlevel = taxlevel, colorPalette = colorIN2018, plotx = "State", facet = "",
             fill_legend = taxlevel, plot_title = "", ylabs = "", xLabs = "State"))

#

(A1 = ggarrange(es2017, lf2017,labels = c("A"), common.legend = T, legend = "right") +
    theme(plot.margin = margin(r=10, l=5)))

(A2 = ggarrange( es2018, lf2018, labels = c("B"), common.legend = T, legend = "right") +
    theme(plot.margin = margin(r=10, l=5)))

ggarrange(A1, A2, ncol = 1)

# --------------------------------------------------------------------------
melt %>%
  dplyr::filter(Year == "2017" & Domain == "Bacteria") %>%
  group_by(Season, State) %>%
  dplyr::mutate(total_abund = sum(Abundance), rel_abund = 100*Abundance / total_abund) %>%
  select(Season, State, !!sym(taxlevel), Abundance, total_abund, rel_abund) %>%
  group_by(Season, State, !!sym(taxlevel)) %>%
  summarise(total = sum(rel_abund)) %>%
  slice_max(order_by = total, n=5) %>%
  write.csv("./output/top5_order_bac_2017_update.csv", quote = F,row.names = F)

melt %>%
  dplyr::filter(Year == "2018" & Domain == "Bacteria") %>%
  group_by(Season, State) %>%
  dplyr::mutate(total_abund = sum(Abundance), rel_abund = 100*Abundance / total_abund) %>%
  select(Season, State, !!sym(taxlevel), Abundance, total_abund, rel_abund) %>%
  group_by(Season, State, !!sym(taxlevel)) %>%
  summarise(total = sum(rel_abund)) %>%
  slice_max(order_by = total, n=5) %>%
  write.csv("./output/top5_order_bac_2018_update.csv", quote = F,row.names = F)


#########################################################################
# BCA 
#########################################################################
rm(list = ls())


# load bca
bca = read.delim("./clean_BCA_bacterial_list.tsv")
head(bca)
bca.sp = unique(bca %>% filter(Functions=="BCA") %>% pull(Species))


#psq2 = readRDS("./bst_psq2_bac.rds")
#phylo.prune.cl.sp = tax_glom(psq2, "Species")
#sp_melt = phylo.prune.cl.sp %>%
#  psmelt()
#head(sp_melt)
#saveRDS(sp_melt, "./sp_melt.rds")



sp_melt = readRDS("./sp_melt.rds")
head(sp_melt)

sp_melt$Species2 = gsub("_", " ", sp_melt$Species)
sp_melt$Species2 = gsub("\\[|\\]", "", sp_melt$Species2)

sp_melt$Species2

(sp_relcomp_season = sp_melt %>%
    dplyr::group_by(Year, State, Season) %>%
    dplyr::mutate(total_abund = sum(Abundance)) %>%
    dplyr::select(Season, State, Species2, Abundance, total_abund) %>%
    dplyr::mutate(rel_abund = 100*Abundance / total_abund)
)

sp_antg = sp_relcomp_season %>%
  dplyr::filter(Species2 %in% bca.sp)

sp_antg


# test on the abundance changes of the most abundance order #--------------------------

# sp_antg %>%
#   group_by(Year, Season, State, Species2) %>%
#   summarise(acumulative_relative_abundance = sum(rel_abund), mean_abundance = mean(rel_abund)) %>%
#   arrange(desc(acumulative_relative_abundance)) %>%
#   write.table("./output/ver4/BCA_relatvie_abundance.tsv", sep="\t", quote = F, row.names = F)

(bca.top10.sp = sp_antg %>%
    group_by(Year, Season, State, Species2) %>%
    summarise(bysum = sum(rel_abund)) %>%
    slice_max(order_by = bysum, n = 10))

head(bca.top10.sp)
bca.top10.sp

bca.top6.unique = unique(bca.top10.sp$Species2)[1:6] # top 6
bca.top6.unique

source("./soilmicrobiome_helpers.R")

taxlevel = "Species2"

bca.kw2017 = kw_test(sp_melt, bac_vector = bca.top6.unique, out = F, outdir = "./output/ver4/bca_kw_test/", taxlevel = taxlevel)

(bca.figurelist2017 = plot_kw(bca.kw2017, bca.top6.unique, ylab = "", taxlevel = taxlevel, facet_ = "~ Species2", gap = 0.01, annoX = 0.5, xlab = ""))

bca.fig.top.2017 = ggpubr::ggarrange(plotlist = bca.figurelist2017$figlist, align = "hv", common.legend = T, legend = "right")

(bca2017 = annotate_figure(bca.fig.top.2017, bottom = text_grob("State"),left = text_grob("Relative abundance (%)", rot = 90, size=12)))

bca2017


#--------------------------------------------------------------------------------------------
bca.kw2018 = kw_test(sp_melt, bac_vector = bca.top6.unique, out = F, outdir = "./output/ver4/bca_kw_test/", taxlevel = taxlevel, filterYear = 2018)

(bca.figurelist2018 = plot_kw(bca.kw2018, bca.top6.unique, ylab = "", taxlevel = taxlevel, facet_ = "~ Species2", gap = 0.03, annoX = 0.5, xlab = ""))

bca.fig.top.2018 = ggarrange(plotlist = bca.figurelist2018$figlist, align = "hv", common.legend = T, legend = "right") 

(bca2018 = annotate_figure(bca.fig.top.2018, bottom = text_grob("State"),left = text_grob("Relative abundance (%)", rot = 90, size=12)))
bca2018


ggarrange(bca2017, bca2018, labels = c("A", "B"), align = "hv", nrow = 2)
