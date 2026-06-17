##############################################################################
##############################################################################
##############################################################################
################ Figures


library(tidyverse) # llamar a la libreria
library(dplyr)
library(perm)
library(devtools)
library(UpSetR)
library("visdat")
library(forcats)
library(ggplot2)
library(ggsignif)
library(tidyr)
library(cowplot)
library(writexl)
library(broom)
library(knitr)    
library(kableExtra)
library(stringr)
library(tibble)
library(cowplot)
library(plotly)

df1 <- read_csv("bd_plt1.csv" )

gene_colors <- c(
  "CD36"   = "#374E55",   # gris pizarra
  "DJ-1"   = "#DF8F44",   # ocre apagado
  "EPG5"   = "#00A1D5",   # azul acero
  "FBXO7"  = "#B24745",   # rojo ladrillo
  "GBA1"   = "#79AF97",   # verde salvia
  "PINK1"  = "#6A6599",   # violeta grisáceo
  "PLA2G6" = "#80796B",   # marrón topo
  "POLG"   = "#D5A823",   # mostaza
  "PRKN"   = "#2D6A9F",   # azul marino
  "SYNJ1"  = "#C08497",   # rosa antiguo
  "VPS13C" = "#3C7A40",   # verde bosque
  "FKBP4"  = "#A8603A",   # marrón cobre
  "LRRK2"  = "#4A4A4A",   # gris carbón
  "PARK7"  = "#8A9A5B",   # verde oliva
  "RET"    = "#5B7FA6"    # azul periwinkle
)



df1 <- df1 %>%
  mutate(
    consequence_group = case_when(
      
      # missing
      consequence == "." ~ "Spl",
      
      # synonymous
      str_detect(consequence, regex("synonymous", ignore_case = TRUE)) &
        !str_detect(consequence, regex("nonsynonymous", ignore_case = TRUE)) ~ "Syn",
      
      # missense
      str_detect(consequence, regex("nonsynonymous", ignore_case = TRUE)) ~ "Mis",
      
      # loss of function
      str_detect(consequence,
                 regex("stopgain|frameshift|startloss",
                       ignore_case = TRUE)) ~ "LoF",
      
      # deletions
      str_detect(consequence, regex("del", ignore_case = TRUE)) ~ "Del",
      
      # duplications
      str_detect(consequence, regex("dup", ignore_case = TRUE)) ~ "Dup",
      
      # inversions
      str_detect(consequence, regex("inv", ignore_case = TRUE)) ~ "Inv",
      
      # fallback
      TRUE ~ "Spl"
    )
  )



df1$consequence_group <- factor(
  df1$consequence_group,
  levels = rev(c(
    "Del",
    "Dup",
    "Inv",
    "Spl",
    "LoF",
    "Syn",
    "Mis"
    
  ))
)



df_sun <- df1 %>%
  mutate(consequence_group = as.character(consequence_group)) %>%
  count(group, gene, type, consequence_group)

# ── Función para construir plot_df por tipo ───────────────────────────────────
make_sunburst <- function(data, tipo) {
  
  df_tipo <- data %>% filter(type == tipo)
  
  # LEVEL 1: group
  lvl1 <- df_tipo %>%
    count(group, wt = n) %>%
    transmute(ids = group, labels = group, parents = "", values = n,
              is_small = FALSE)
  
  # LEVEL 2: gene (hijo de group)
  lvl2 <- df_tipo %>%
    count(group, gene, wt = n) %>%
    transmute(
      ids      = paste(group, gene, sep = "-"),
      labels   = gene,
      parents  = group,
      values   = n,
      is_small = FALSE
    )
  
  # LEVEL 3: consequence_group + N
  lvl3 <- df_tipo %>%
    count(group, gene, consequence_group, wt = n) %>%
    filter(!is.na(consequence_group),
           trimws(as.character(consequence_group)) != "") %>%
    mutate(
      pct      = n / sum(n) * 100,
      is_small = pct < 2,
      labels   = ifelse(is_small,
                        as.character(consequence_group),
                        paste0(consequence_group, " (", n, ")"))
    ) %>%
    transmute(
      ids      = paste(group, gene, consequence_group, sep = "-"),
      labels   = labels,
      parents  = paste(group, gene, sep = "-"),
      values   = n,
      is_small = is_small
    )
  
  # Unir y colorear
  plot_df <- bind_rows(lvl1, lvl2, lvl3) %>%
    mutate(
      gene_color = case_when(
        grepl("CD36",   ids) ~ "#374E55",
        grepl("DJ-1",   ids) ~ "#DF8F44",
        grepl("EPG5",   ids) ~ "#00A1D5",
        grepl("FBXO7",  ids) ~ "#B24745",
        grepl("GBA1",   ids) ~ "#79AF97",
        grepl("PINK1",  ids) ~ "#6A6599",
        grepl("PLA2G6", ids) ~ "#80796B",
        grepl("POLG",   ids) ~ "#D5A823",
        grepl("PRKN",   ids) ~ "#2D6A9F",
        grepl("SYNJ1",  ids) ~ "#C08497",
        grepl("VPS13C" , ids) ~ "#3C7A40",
        grepl("FKBP4",  ids) ~ "#A8603A",
        grepl("LRRK2",  ids) ~ "#4A4A4A",
        grepl("PARK7",  ids) ~ "#8A9A5B",
        grepl("RET",    ids) ~ "#5B7FA6",
        TRUE ~ "#D3D3D3"
      )
    )
  
  return(plot_df)
}

# ── Función para hacer el plot_ly ─────────────────────────────────────────────
make_plot <- function(plot_df, titulo) {
  
  colors_vec <- case_when(
    plot_df$ids == "CH" ~ "#000000",   # dark blue
    plot_df$ids == "DI" ~ "#4b4b4b",   # dark red
    TRUE ~ plot_df$gene_color
  )
  
  plot_df <- plot_df %>%
    mutate(
      label_plot = ifelse(is_small, "", labels)
    )
  
  plot_ly(
    data          = plot_df,
    ids           = ~ids,
    labels        = ~label_plot,
    parents       = ~parents,
    values        = ~values,
    type          = "sunburst",
    branchvalues  = "total",
    textinfo      = NA,
    hovertemplate = ~paste0(
      labels, " (", values, ")",
      "<extra></extra>"
    ),
    textfont = list(color = "white", size = 16),
    marker = list(
      colors = colors_vec,
      line   = list(color = "white", width = 1)
    )
  ) %>%
    layout(
      title  = list(text = paste0("<b>", titulo, "</b>"),
                    font = list(size = 18)),
      margin = list(b = 30, l = 30, r = 30, t = 60)
    )
}

# ── Construir los 2 plot_df ───────────────────────────────────────────────────
df_snvs <- make_sunburst(df_sun, "SNVs")
df_sv   <- make_sunburst(df_sun, "SV")

# ── Plotear ───────────────────────────────────────────────────────────────────
p_snvs <- make_plot(df_snvs, "SNVs")
p_sv   <- make_plot(df_sv,   "SV")

plot_df <- df_snvs %>%
  mutate(
    is_outer = str_count(ids, "-") == 2,
    label_plot = ifelse(is_outer, "", labels)
  )

p_snvs
p_sv


library(htmlwidgets)
library(webshot2)

saveWidget(p_snvs, "snv_tmp.html", selfcontained = TRUE)
saveWidget(p_sv, "sv_tmp.html", selfcontained = TRUE)

webshot(
  "snv_tmp.html",
  "SNV_sunburst.png",
  vwidth = 1000,
  vheight = 1000,
  zoom = 3
)

webshot(
  "sv_tmp.html",
  "SV_sunburst.png",
  vwidth = 1000,
  vheight = 1000,
  zoom = 3
)


#########################################
####### Figure 2A

df1 <- read_csv("method_ancestry.csv" )


df1 %>%
  count(method, sort = TRUE) %>%
  print(n = Inf)   # muestra todos sin truncar

df1 <- df1 %>%
  mutate(
    method_group = case_when(
      
      # ── MIX: más de un método (va PRIMERO) ──────────────────────────────────
      grepl("[,;]|\\band\\b|\\+", method, ignore.case = TRUE)
      ~ "Mix",
      
      # ── MLPA solo ────────────────────────────────────────────────────────────
      grepl("^MLPA$", method, ignore.case = TRUE)
      ~ "MLPA",
      
      # ── WGS ──────────────────────────────────────────────────────────────────
      grepl("WGS", method, ignore.case = TRUE) |
        grepl("paired-end sequencing", method, ignore.case = TRUE)
      ~ "WGS",
      
      # ── WES ──────────────────────────────────────────────────────────────────
      grepl("WES", method, ignore.case = TRUE)
      ~ "WES",
      
      # ── Panel NGS ────────────────────────────────────────────────────────────
      grepl("panel", method, ignore.case = TRUE) |
        grepl("targeted sequenc", method, ignore.case = TRUE) |
        grepl("targetted sequenc", method, ignore.case = TRUE) |
        grepl("Multigene", method, ignore.case = TRUE) |
        grepl("hybrid-capture", method, ignore.case = TRUE) |
        grepl("High-throughout sequencing", method, ignore.case = TRUE)
      ~ "Panel NGS",
      
      # ── Sanger ───────────────────────────────────────────────────────────────
      grepl("^Sanger", method, ignore.case = TRUE) |
        grepl("^direct sequencing$", method, ignore.case = TRUE) |
        grepl("^sequencing$", method, ignore.case = TRUE) |
        grepl("^PCR-amplified and directly sequenced", method, ignore.case = TRUE)
      ~ "Sanger",
      
      # ── PCR-based ────────────────────────────────────────────────────────────
      grepl("PCR", method, ignore.case = TRUE) |
        grepl("DHPLC", method, ignore.case = TRUE) |
        grepl("RT-PCR", method, ignore.case = TRUE) |
        grepl("RFLP", method, ignore.case = TRUE) |
        grepl("TaqMan", method, ignore.case = TRUE) |
        grepl("SSCP", method, ignore.case = TRUE) |
        grepl("STR genotyping", method, ignore.case = TRUE) |
        grepl("Transgenomic", method, ignore.case = TRUE)
      ~ "PCR-based",
      
      # ── Funcional / Otro ─────────────────────────────────────────────────────
      grepl("CRISPR|iPSC|Western|Inmunof|ELISA|Flow Cyto|RNA analysis|Inmunoblot|Sendai",
            method, ignore.case = TRUE)
      ~ "Funcional/Otro",
      
      # ── No reportado ─────────────────────────────────────────────────────────
      is.na(method)
      ~ "Not Reported",
      
      TRUE ~ "Otro"
    )
  )

# ── Verificar ─────────────────────────────────────────────────────────────────
df1 %>%
  count(method_group, method, sort = TRUE) %>%
  print(n = Inf)

df1 <- df1 %>%
  mutate(
    ancestry_group = case_when(
      
      # EAS - East Asian
      ancestry %in% c("Chinese", "Taiwanese (ethnic Chinese)", "China", 
                      "Japanese", "Korean", "Vietnam", "Taiwan", 
                      "South Korea", "Japan", "Filipino")
      ~ "EAS",
      
      # EUR - European
      ancestry %in% c("Belgian", "Caucasian", "Italian", "Irish", "Czech",
                      "Polish", "Ukrainian", "Swiss", "Netherlands",
                      "European", "Serbian", "Norway", "Spain", "Italy",
                      "French", "France", "UK", "Finnland", "North American")
      ~ "EUR",
      
      # AFR - African
      ancestry %in% c("African", "South African", "Tunisia")
      ~ "AFR",
      
      # AMR - American
      ancestry %in% c("Mexican", "Brazil", "Mexico", "Peru")
      ~ "AMR",
      
      # MDE - Middle East
      ancestry %in% c("Yemen", "Turkish", "Cyprus")
      ~ "MDE",
      
      is.na(ancestry) ~ "No reportado",
      TRUE ~ "Otro"
    )
  )

# ── Verificar ─────────────────────────────────────────────────────────────────
df1 %>% count(ancestry_group, ancestry, sort = TRUE) %>% print(n = Inf)

df1 %>%
  count(ancestry_group) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  arrange(desc(n))


# ── Paleta por method_group ───────────────────────────────────────────────────
method_colors <- c(
  "WGS"            = "#1B3A6BFF",   # azul noche
  "WES"            = "#4E79A5FF",   # azul acero
  "Panel NGS"      = "#76A1CDFF",   # azul claro
  "Sanger"         = "#FFDE90FF",   # dorado arena
  "MLPA"           = "#F8CA7CFF",   # café arena
  "Mix"    = "#CE8A37FF",   # café topo
  "PCR-based"      = "#AC420AFF",   # gris medio
  "Not Reported" = "#752305FF"  # beige  
)

# ── Preparar datos ────────────────────────────────────────────────────────────
df_stack <- df1 %>%
  filter(!is.na(ancestry_group), ancestry_group != "No reportado") %>%
  count(ancestry_group, method_group) %>%
  group_by(ancestry_group) %>%
  mutate(pct = round(n / sum(n) * 100, 1)) %>%
  ungroup() %>%
  mutate(
    ancestry_group = factor(ancestry_group, levels = c("AFR","AMR","EAS","EUR","MDE")),
    method_group   = factor(method_group, levels = names(method_colors))
  )

# ── Plot ──────────────────────────────────────────────────────────────────────
p3 <- ggplot(df_stack,
             aes(x    = ancestry_group,
                 y    = pct,
                 fill = method_group)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(pct >= 5, paste0(pct, "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 5, color = "white", fontface = "bold") +
  scale_fill_manual(values = method_colors,
                    name   = "Method") +
  scale_y_continuous(expand = c(0, 0),
                     labels = scales::percent_format(scale = 1)) +
  labs(
    x        = NULL,
    y        = "% variants"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position    = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(face = "bold", size = 12),
    axis.text.x        = element_text(color = "gray50", size = 12)
  )
ggsave("Ancestry_stacked.png", plot = p3, width = 9, height = 7, dpi = 300)
p3 


######################## Figure 2B

library(ggalluvial)
library(ggbreak)

df <- read_csv("variant_com_df.csv" )

df <- df %>%
  rename(
    variant1 =`variant 1`,
    variant2 = `variant 2`
  )


df$variant1_type <- case_when(
  str_detect(df$variant1, regex("Ex", ignore_case = TRUE)) ~ "SV",
  str_detect(df$variant1, regex("p\\.|c\\.", ignore_case = TRUE)) ~ "SNV",
  TRUE ~ NA_character_
)


df$variant2_type <- case_when(
  str_detect(df$variant2, regex("Ex", ignore_case = TRUE)) ~ "SV",
  str_detect(df$variant2, regex("p\\.|c\\.", ignore_case = TRUE)) ~ "SNV",
  TRUE ~ NA_character_
)

df <- df %>%
  mutate(
    combination = case_when(
      
      # SNV + SNV
      variant1_type == "SNV" & variant2_type == "SNV" ~ "SNV+SNV",
      
      # SV + SV
      variant1_type == "SV" & variant2_type == "SV" ~ "SV+SV",
      
      # combinaciones mixtas
      (variant1_type == "SNV" & variant2_type == "SV") |
        (variant1_type == "SV" & variant2_type == "SNV") ~ "SNV+SV",
      
      TRUE ~ NA_character_
    )
  )

ggplot(df,
       aes(x = Gene,
           fill = combination)) +
  
  geom_bar(position = "stack",
           color = "black",
           width = 0.8) +
  
  scale_y_continuous(labels = scales::percent) +
  
  scale_fill_manual(values = c(
    "SNV+SNV" = "#4E79A7",
    "SNV+SV"  = "#F28E2B",
    "SV+SV"   = "#8C564B"
  )) +
  
  geom_text(stat = "count",
            aes(label = after_stat(count)),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white", fontface = "bold")
labs(
  x = "Gene",
  y = "Percentage of subjects",
  fill = "Variant combination"
) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      face = "bold"
    ),
    
    legend.position = "top"
  )

ggplot(df,
       aes(x = Gene,
           fill = combination)) +
  geom_bar(position = "stack",
           color = "black",
           width = 0.8) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c(
    "SNV+SNV" = "#4E79A7",
    "SNV+SV"  = "#F28E2B",
    "SV+SV"   = "#8C564B"
  )) +
  geom_text(stat = "count",
            aes(label = after_stat(count)),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white", fontface = "bold") +
  labs(
    x = "Gene",
    y = "Number of subjects",
    fill = "Variant combination"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),   # <- quita cuadrícula horizontal
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position    = "top"
  )

# ── Agrupar genes: PRKN, PINK1 y Otros ───────────────────────────────────────
df_plot <- df %>%
  filter(!is.na(combination)) %>%
  mutate(
    gene_group = case_when(
      Gene == "PRKN"  ~ "PRKN",
      Gene == "PINK1" ~ "PINK1",
      TRUE            ~ "Others"
    ),
    gene_group = factor(gene_group, levels = c("PRKN", "PINK1", "Others"))
  ) %>%
  count(gene_group, combination) %>%
  group_by(gene_group) %>%
  mutate(n_subjects = n) %>%
  ungroup()

# ── Plot ──────────────────────────────────────────────────────────────────────
p2 <- ggplot(df_plot,
             aes(x    = gene_group,
                 y    = n_subjects,
                 fill = combination)) +
  geom_col(position = "stack",
           color = "white",
           linewidth = 0.4,
           width = 0.7) +
  scale_fill_manual(values = c(
    "SNV+SNV" = "#78518b",
    "SNV+SV"  = "#265596",
    "SV+SV"   = "#b81f30"
  )) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x    = "Gene",
    y    = "Number of subjects",
    fill = "Variant combinations"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(face = "bold", size = 12),
    axis.text.y        = element_text(color = "gray50", size = 12),
    legend.position    = "top"
  ) +
  scale_y_break(c(45, 160))

ggsave("Variant_comb2.png", plot = p2, width = 9, height = 7, dpi = 300)
p2

genes_others <- df %>%
  filter(!is.na(combination)) %>%
  filter(!Gene %in% c("PRKN", "PINK1")) %>%
  distinct(Gene) %>%
  arrange(Gene) %>%
  pull(Gene)

genes_others
