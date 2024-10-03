
  

library(dplyr)
library(phylotools)
library(stringr)
library(tidyr)
library(ggplot2)
library(readxl)
rm(list = ls())

# Manipulate fasta alignment of sequences.

# See if other organisms match HS sequence

setwd("C:/users/ezraa/desktop/plots 9.5.24/text_plots/")
fasta = phylotools::read.fasta("C:/Users/ezraa/Desktop/plots 9.5.24/alignment/pnpo_cobalt_alignment.fa") %>% 
  mutate(name = str_extract(seq.name, pattern = '\\[organism=.*\\]') %>% gsub(pattern = "\\].*|\\[organism=", replacement = "")) %>%
  # mutate(len = nchar(seq.text)) %>% 
  select(-seq.name)


fastaw = fasta %>% select(-name) %>% t() %>% as.data.frame() %>% setNames(fasta$name) %>% mutate(across(everything(), ~ str_split(.,pattern = ""))) %>% unnest(cols = everything())

fastaw =fastaw %>%  
  mutate(whole_index = 1:264) %>%
  mutate(across(-whole_index, 
                list(index = ~ if_else(.!="-", 
                                       cumsum(.!="-"), # Cumulative sum to create index
                                       NA_integer_)), 
                .names = "{.col}_index")) 

colnames(fastaw)[1:6] = paste0(colnames(fastaw)[1:6], "_aa")



aa = fastaw %>% select(ends_with("_aa"),whole_index) %>% pivot_longer(-whole_index, values_to = "aa") %>% mutate(name = gsub(name, pattern = "_aa", replacement = ""))

index = fastaw %>% select(ends_with("_index"),whole_index) %>% pivot_longer(-whole_index, values_to = "index") %>% mutate(name = gsub(name, pattern = "_index", replacement = ""))

fastal = merge(aa, index)


organisms = fastaw %>% select(-whole_index) %>% colnames() %>% gsub(pattern = "_.*", replacement = "") %>% unique()
paste0(organisms, sep = "", collapse = "|")

rm(fasta, fastaw, aa, index)

# Now do binding sites


fmn = read_xlsx("C:/Users/ezraa/Desktop/plots 9.5.24/musayev binding sites.xlsx", sheet = "fmn") %>% mutate("ligand" = "FMN")
plp = read_xlsx("C:/Users/ezraa/Desktop/plots 9.5.24/musayev binding sites.xlsx", sheet = "plp") %>% mutate("ligand" = "PLP")

ligands = rbind(fmn, plp) %>% rename(ligand_atom = Ligand, protein = Proteinb)
rm(fmn, plp)

ligands = ligands %>% mutate(index = gsub(protein, pattern = "[A-Za-z]{3}|[A-Za-z]{1}\\s.*", replacement = "") %>% as.numeric()) %>% select(-c(ligand_atom, Water, protein)) %>% pivot_longer(cols = -c(ligand,index)) %>% filter(value != "	
—") %>% select(-value) %>% unique()

# combine ligand and aa information


fastaf = left_join(fastal, ligands) 
rm(fastal)


fastaf$match = NA
for(i in 1:nrow(fastaf)){
  onerow = fastaf[i,]
  
  if(onerow$name == "Homo sapiens"){
    fastaf$match[i] = "Match"
  } else if(onerow$aa == fastaf$aa[fastaf$name=="Homo sapiens" & fastaf$whole_index == onerow$whole_index]){
    fastaf$match[i] = "Match"
  } else{
    fastaf$match[i] = "Mismatch"
  }
}
rm(onerow,ind,i)

# Now make plots


my_theme = theme(
  panel.grid = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.y = element_text(size = 15),
  axis.text.x = element_text(size = 15),
  # legend.position = 'none',
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0.5, size = 40)
) 

match_color = "black"
mismatch_color = "grey"

fmn_color = "#00A08A"
plp_color = "#FC1A00"

label_index = fastaf$index[order(fastaf$whole_index)][fastaf$name == "Homo sapiens"]

# Mutations


mutations = read_xlsx("C:/Users/ezraa/Desktop/PNPO 7.9.24/final PNPO.xlsx") %>%
  # select only homozygous missense
  filter(grepl(Type, pattern = "Missense|missense")==TRUE)


table(mutations$Type)


hom = mutations[mutations$`2nd protein mutation`=="hom",] %>% 
  select(`1st protein mutation`) %>%
  rename(mutation = `1st protein mutation`)
het = mutations[mutations$`2nd protein mutation`!="hom",] %>% 
  select(c(`1st protein mutation`, `2nd protein mutation`)) %>% 
  pivot_longer(cols = everything(), values_to = "mutation") %>% 
  filter(grepl("deletion|fs|c.|del", mutation)==FALSE) %>%
  separate_longer_delim(cols = mutation, delim = "+") %>%
  select("mutation") %>%
  mutate(mutation = trimws(mutation))

imdonefunction = function(x){
  if(length(x) > 1){
    return(1)
  }
else{
  return(0)
}
}
mutations = rbind(het, hom) %>% 
  mutate(index = gsub("[A-Za-z]", replacement = "",mutation) %>% str_split("_")) 

# mutations$index2 = sapply(mutations$index, FUN = imdonefunction)



rm(het,hom)

mutations$index = as.numeric(mutations$index)
mutations$name = ""

mutations = unique(mutations)
index_df = fastaf[fastaf$name=="Homo sapiens",c("whole_index", "index")]

mutations = left_join(mutations, index_df)
rm(index_df)

mutations$aa = "*"
mutations$match = "Match"
mutations$name = " "


fastaf = full_join(fastaf, mutations)
fastaf$name = factor(fastaf$name, levels = c("Rattus norvegicus","Mus musculus", "Bos taurus","Saccharomyces cerevisiae S288C","Escherichia coli K-12","Homo sapiens", " "), ordered = TRUE)


unique(fastaf$name)

# Plots:
library(cowplot) 

y_label_plot = ggplot() + 
  geom_text(data = fastaf, size = 15, aes(x = whole_index, y = name, label = aa, color = match)) +
  geom_tile(data = fastaf, aes(x = whole_index, y = name, fill = ligand)) +
  my_theme + 
  theme(legend.position = "left", axis.text.y = element_text(face = "bold", color = "black", size = 20)) +
  scale_y_discrete(labels = function(y) str_wrap(y, width = 30)) + 
  scale_color_manual(name = "Human PNPO alignment match", values = c(match_color, mismatch_color)) +
  scale_fill_manual(name = "Ligand", na.translate = FALSE, na.value = "transparent", values = c(plp_color, fmn_color)) +
  xlim(5,10)


y_labels = get_plot_component(y_label_plot, "axis-l", "ylab") 
y_labels = ggdraw(y_labels)

legend = get_plot_component(y_label_plot, "guide-box-left")
legend = ggdraw(legend)

legend
ggsave(legend, filename = "legend.tif", height = 7)

y_labels
bat = theme(axis.text.y = element_blank(), 
                         axis.ticks.y = element_blank(),
                         axis.title.y = element_blank(),
                         legend.position = "none",
                         panel.border = element_rect(color = "black", fill = NA, size = 2))

# get rid of NAs
label_index[is.na(label_index)] = ""

  # Asp33Val

asp33 = ggplot() + 
  geom_text(data = fastaf, size = 15, aes(x = whole_index, y = name, label = aa, color = match)) +
  geom_tile(data = fastaf, aes(x = whole_index, y = name, fill = ligand)) +
  scale_x_continuous(breaks = seq(31, 35, by = 1), 
                     limits = c(31,35), labels = label_index[31:35], 
                     expand = expansion(add = 1)) +
  # coord_fixed(ratio = 1) +
  scale_color_manual(values = c(match_color, mismatch_color)) +
  scale_fill_manual(na.value = "transparent", values = c(plp_color, fmn_color)) + 
  my_theme 

# ggsave(asp33, filename = "asp33.png")

glu50leu53 = ggplot() + 
  geom_text(data = fastaf, size = 15, aes(x = whole_index, y = name, label = aa, color = match)) +
  geom_tile(data = fastaf, aes(x = whole_index, y = name, fill = ligand)) +
  scale_x_continuous(breaks = seq(47, 56, by = 1), 
                     limits = c(47,56), labels = label_index[47:56], 
                     expand = expansion(add = 1)) +
  # coord_fixed(ratio = 1) +
  scale_color_manual(values = c(match_color, mismatch_color)) +
  scale_fill_manual(na.value = "transparent", values = c(plp_color, fmn_color)) + 
  my_theme 

# ggsave(glu50leu53, filename = "glu50leu53.png")



asp33glu50leu53 = plot_grid(y_labels,
                            asp33+bat, 
                            glu50leu53+bat, 
                            nrow = 1, 
                            rel_widths = c(5,7,10),
                            align = "h")

asp33glu50leu53

# ggsave(asp33glu50leu53, filename = "asp33glu50leu53.png")
# ala78


ala78 = ggplot() + 
  geom_tile(data = fastaf, alpha = .5, aes(x = whole_index, y = name, fill = ligand)) +
  geom_text(data = fastaf, size = 15, aes(x = whole_index, y = name, label = aa, color = match)) +
  scale_x_continuous(breaks = seq(76, 80, by = 1), 
                     limits = c(76,80), 
                     labels = label_index[76 :80],
                     expand = expansion(add = 1)) +
  # coord_fixed(ratio = .8) + 
  scale_color_manual(values = c(match_color, mismatch_color)) +
  scale_fill_manual(na.value = "transparent", values = c(plp_color, fmn_color)) +
  my_theme 

# ggsave(ala78, filename = "ala78.png")

# arg95


arg95 = ggplot() + 
  geom_tile(data = fastaf, alpha = .5, aes(x = whole_index, y = name, fill = ligand)) +
  geom_text(data = fastaf, size = 15, aes(x = whole_index, y = name, label = aa, color = match)) +
  scale_x_continuous(breaks = seq(94, 98, by = 1), 
                     limits = c(94,98), 
                     labels = label_index[94:98],  expand = expansion(add = 1)) +
  # coord_fixed(ratio = .8) + 
  scale_color_manual(values = c(match_color, mismatch_color)) +
  scale_fill_manual(na.value = "transparent", values = c(plp_color, fmn_color)) +
  my_theme 

# ggsave(arg95, filename = "arg95.png")
# plp_binding_site


plp_binding_site = ggplot() + 
  geom_tile(data = fastaf, alpha = .5, aes(x = whole_index, y = name, fill = ligand)) +
  geom_text(data = fastaf, size = 15, aes(x = whole_index, y = name, label = aa, color = match)) +
  scale_x_continuous(breaks = seq(115, 129, by = 1), 
                     limits = c(115,129), 
                     labels = label_index[115:129], 
                     expand = expansion(add = 1)) +
  # coord_fixed(ratio = ) +
  scale_color_manual(values = c(match_color, mismatch_color)) +
  scale_fill_manual(na.value = "transparent", values = c(plp_color, fmn_color)) +
  my_theme 


plp_binding_site

# ggsave(plp_binding_site, filename = "plp_binding_site.png")

# arg138arg141


arg138arg141 = ggplot() + 
  geom_tile(data = fastaf, alpha = .5, aes(x = whole_index, y = name, fill = ligand)) +
  geom_text(data = fastaf, size = 15, aes(x = whole_index, y = name, label = aa, color = match)) +
  scale_x_continuous(breaks = seq(138, 147, by = 1), 
                     limits = c(138,147), 
                     labels = label_index[138:147], 
                     expand = expansion (add = 1)) +
  # coord_fixed(ratio = .8) + 
  scale_color_manual(values = c(match_color, mismatch_color)) +
  scale_fill_manual(na.value = "transparent", values = c(plp_color, fmn_color)) +
  my_theme 

# ggsave(arg138arg141, filename = "arg138arg141.png")
# arg161ile167


arg161ile167 = ggplot() + 
  geom_tile(data = fastaf, alpha = .5, aes(x = whole_index, y = name, fill = ligand)) +
  geom_text(data = fastaf, size = 15, aes(x = whole_index, y = name, label = aa, color = match)) +
  scale_x_continuous(breaks = seq(158, 172, by = 1), limits = c(158,172), labels = label_index[158:172], expand = expansion(add = 1)) +
  # coord_fixed(ratio = .8) + 
  scale_color_manual(values = c(match_color, mismatch_color)) +
  scale_fill_manual(na.value = "transparent", values = c(plp_color, fmn_color)) +
  my_theme 

# ggsave(arg161ile167,filename = "arg161ile167.png")
# arg161

# pro213


pro213 = ggplot() + 
  geom_tile(data = fastaf, alpha = .5, aes(x = whole_index, y = name, fill = ligand)) +
  geom_text(data = fastaf, size = 15, aes(x = whole_index, y = name, label = aa, color = match)) +
  scale_x_continuous(breaks = seq(214, 218, by = 1), 
                     limits = c(214,218), 
                     labels = label_index[214:218], 
                     expand = expansion(add = 1)) +
  # coord_fixed(ratio = .8) + 
  scale_color_manual(values = c(match_color, mismatch_color)) +
  scale_fill_manual(na.value = "transparent", values = c(plp_color, fmn_color)) +
  my_theme 


# ggsave(pro213, filename = "pro213.png")

# arg225arg229


arg225arg229 = ggplot() + 
  geom_tile(data = fastaf, alpha = .5, aes(x = whole_index, y = name, fill = ligand)) +
  geom_text(data = fastaf, size = 15, aes(x = whole_index, y = name, label = aa, color = match)) +
  scale_x_continuous(breaks = seq(226, 234, by = 1), 
                     limits = c(226,234), 
                     labels = label_index[226:234]) +
  coord_fixed(ratio = .8) + 
  scale_color_manual(values = c(match_color, mismatch_color)) +
  scale_fill_manual(na.value = "transparent", values = c(plp_color, fmn_color)) +
  my_theme 

# arg225arg229arg234


arg225arg229arg234 = ggplot() + 
  geom_tile(data = fastaf, alpha = .5, aes(x = whole_index, y = name, fill = ligand)) +
  geom_text(data = fastaf, size = 15, aes(x = whole_index, y = name, label = aa, color = match)) +
  scale_x_continuous(breaks = seq(226, 240, by = 1), 
                     limits = c(226,240), 
                     labels = label_index[226:240],
                     expand = expansion(add = 1)) +
  # coord_fixed(ratio = .8) + 
  scale_color_manual(values = c(match_color, mismatch_color)) +
  scale_fill_manual(na.value = "transparent", values = c(plp_color, fmn_color)) +
  my_theme 

ggsave(arg225arg229arg234, filename = "arg225arg229arg234.png")

# tyr254


tyr254 = ggplot() + 
  geom_tile(data = fastaf, alpha = .5, aes(x = whole_index, y = name, fill = ligand)) +
  geom_text(data = fastaf, size = 15, aes(x = whole_index, y = name, label = aa, color = match)) +
  scale_x_continuous(breaks = seq(255, 259, by = 1), 
                     limits = c(255,259), 
                     labels = label_index[255:259],
                     expand = expansion(add = 1)) +
  # coord_fixed(ratio = .8) + 
  scale_color_manual(values = c(match_color, mismatch_color)) +
  scale_fill_manual(na.value = "transparent", values = c(plp_color, fmn_color)) +
  my_theme 

# windows()
tyr254
# ggsave(tyr254, filename = "tyr254.png")


# windows()
# asp33glu50leu53 = plot_grid(y_labels, asp33+bat, glu50leu53+bat, nrow = 1, rel_widths = c(1, 1, 1*(11/6)), rel_heights = c(1,1))
# 
# ggsave(asp33glu50leu53, filename = "C:/users/ezraa/desktop/asp33glu50leu53.png", 
#        dpi = 300,
#        width = 20, 
#        height = 10)
# 
# ggsave(plp_binding_site, filename = "C:/users/ezraa/desktop/plp_binding_site.png", 
#        dpi = 300,
#        width = 20, 
#        height = 10)

# windows()
# plot_grid(
#   # y_labels,
#           asp33+bat, 
#           glu50leu53+bat, 
#           ala78+bat,
#           arg95+bat,
#           # y_labels,
#           plp_binding_site+bat,
#           arg138arg141+bat, 
#           arg161ile167+bat,
#           pro213+bat,
#           arg225arg229arg234+bat,
#           tyr254+bat,
#           rel_widths = c(1,2,1,1,3,2,3,1,3,1)
# )


# row1 = plot_grid(
#   y_labels,
#   asp33+bat, 
#   glu50leu53+bat, 
#   ala78+bat,
#   arg95+bat,
#   rel_widths = c(1,1,2,1,1),
#   nrow = 1
# )


row1 = plot_grid(
  asp33+bat, 
  glu50leu53+bat, 
  ala78+bat,
  arg95+bat,
  plp_binding_site+bat,
  arg138arg141+bat,
  rel_widths = c(1,2,1,1,3,2),
  rel_heights = c(1,1,1,1,1,2),
  nrow = 1, 
  align = "v"
)
ggsave(row1, filename = "row1.jpeg", width = 50, height = 7, limitsize = FALSE)

row2 = plot_grid(
          arg161ile167+bat,
          pro213+bat,
          arg225arg229arg234+bat,
          tyr254+bat,
          rel_widths = c(3,1,3,1),
  nrow = 1, 
  align = "v"
)

ggsave(row2, filename = "row2.jpeg", width = 40, height = 7, limitsize = FALSE)

ggsave(y_labels, filename = "ylabels.jpeg", width = 5, height = 7, limitsize = FALSE)

