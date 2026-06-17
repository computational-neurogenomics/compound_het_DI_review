library(dplyr)
library(stringr)
library(stringi)
library(fuzzyjoin)
library(stringr)
library(readr)
library(tidyr)
library(ggplot2)

# ==========================================================
# 1. Literature Search Update (2025 -> 2026)
# ==========================================================
di_raw <- read.csv("digenic_list_march31.csv", sep = ",")



di_raw$var_1_g <- paste0(di_raw$Gene_1, "_", di_raw$var_1_norm)
di_raw$var_2_g <- paste0(di_raw$Gene_2, "_", di_raw$var_2_norm)

di_raw <- di_raw[, !grepl("^X", names(di_raw))]

# ==========================================================
# 2. Load Digenic Cases
# ==========================================================

unique_variants <- unique(c(di_raw$var_1_g, di_raw$var_2_g))
table(di_raw$Gene_1)
table(di_raw$Gene_2)

length(unique_variants)

combo_counts_gene <- di_raw %>%
  rowwise() %>%
  mutate(
    gene_a = sort(c(Gene_1, Gene_2))[1],
    gene_b = sort(c(Gene_1, Gene_2))[2]
  ) %>%
  ungroup() %>%
  count(gene_a, gene_b, sort = TRUE)

combo_counts_gene



# ==========================================================
# 3. Variant Annotation (ClinVar + HGMD + Genoox)
# ==========================================================

###################### Variant annotation USING CLINVAR
is_rearrangement <- grepl("_(Ex|ex).*(del|dup|inv|partial)", unique_variants)

unique_variants <- unique_variants[!is_rearrangement]



vars_refseq <- as.data.frame(unique_variants)

clinvar_dir <- "clinvar_tables"

annotate_variant <- function(var) {
  
  gene <- sub("_.*", "", var)
  query <- sub(".*_[pc]\\.", "", var)
  
  file <- file.path(clinvar_dir, paste0(gene, ".txt"))
  
  if (!file.exists(file)) {
    return(data.frame(
      unique_variants = var,
      gene = gene,
      query = query,
      match_status = "gene_file_not_found",
      stringsAsFactors = FALSE
    ))
  }
  
  df <- read_tsv(
    file,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    name_repair = "unique"
  )
  
  hits <- df[grep(query, df$Name, fixed = TRUE), , drop = FALSE]
  
  if (nrow(hits) == 0) {
    return(data.frame(
      unique_variants = var,
      gene = gene,
      query = query,
      match_status = "no_match",
      stringsAsFactors = FALSE
    ))
  }
  
  hits$unique_variants <- var
  hits$gene <- gene
  hits$query <- query
  hits$match_status <- ifelse(nrow(hits) == 1, "matched", "multiple_matches")
  
  hits
}

results <- bind_rows(lapply(vars_refseq$unique_variants, annotate_variant))



hgmd_df <- read.csv("NBA_VARS.csv") #### Supplementary table 5 from DOI: 10.1101/2023.11.06.23298176


vars_clean <- vars_refseq %>%
  mutate(
    gene = sub("_.*", "", unique_variants),
    aa = sub(".*_p\\.", "", unique_variants),
    
    aa_from = stringr::str_extract(aa, "^[A-Za-z]{3}"),
    aa_to   = stringr::str_extract(aa, "[A-Za-z]{3}$"),
    codon   = as.numeric(stringr::str_extract(aa, "[0-9]+")),
    
    is_trunc = stringr::str_detect(aa, "fs|Ter|\\*|X$"),
    
    aa_change = ifelse(
      is_trunc,
      paste0(aa_from, "-Term"),
      paste0(aa_from, "-", aa_to)
    )
  )

hgmd_clean <- hgmd_df %>%
  rename(
    gene = `Gene.symbol`,
    aa_change = `Amino.Acid.change`,
    codon = `Codon.number`,
    dbsnp = dbsnp,
    hgvs = HGVS
  ) %>%
  mutate(
    codon = as.numeric(codon)
  )

annotated <- vars_clean %>%
  left_join(hgmd_clean, by = c("gene", "aa_change", "codon"))

table(is.na(annotated$hgvs))

round2 <- merge(results, annotated, by = "unique_variants")

round2 <- round2 %>%
  mutate(
    captured = match_status=="matched" | !is.na(hgvs)
  )

table(round2$captured)

round2_noCNV <- round2[!grepl("Ex", round2$unique_variants), ]
table(round2_noCNV$captured)

not_annotated <- round2_noCNV %>%
  filter(match_status!="matched" & is.na(hgvs))


write.csv(not_annotated, 
          "not_annotated_variants_DIGENIC_RAW.csv", 
          row.names = FALSE)



############## Manual annotated variants using Genoox service Franklin: https://franklin.genoox.com/

genoox_raw <- read.csv("not_annotated_variants_DIGENIC.csv")

round3 <- merge(round2_noCNV, genoox_raw, by = "unique_variants", all = TRUE)

toannot <- round3
colnames(toannot)




# ==========================================================
# 4. CPRA Harmonisation
# ==========================================================


################### Parsing columns from different sources to get CPRA 


# ----------------------------
# 1) helper: complement bases
# ----------------------------
comp_base <- function(x) {
  x <- toupper(x)
  chartr("ACGT", "TGCA", x)
}

# ----------------------------
# 2) ClinVar SPDI -> CPRA
# example:
# NC_000006.12:162443313:A:T  ->  6:162443314:A:T
# IMPORTANT: +1 to position
# ----------------------------
spdi_to_cpra <- function(spdi) {
  if (is.na(spdi) || spdi == "") return(NA_character_)
  
  parts <- str_split(spdi, ":", simplify = TRUE)
  if (ncol(parts) < 4) return(NA_character_)
  
  seqid <- parts[1]
  pos0  <- suppressWarnings(as.numeric(parts[2]))
  ref   <- parts[3]
  alt   <- parts[4]
  
  chr <- str_match(seqid, "^NC_0*([0-9XYM]+)\\.[0-9]+$")[,2]
  if (is.na(chr) || is.na(pos0)) return(NA_character_)
  
  pos1 <- pos0 + 1
  
  paste(chr, pos1, ref, alt, sep = ":")
}

# ----------------------------
# 3) NBA coords + Base.change -> CPRA
# examples:
# Genomic.coordinates..GRCh38. = "chr6:161386823 (-)"
# Base.change = "C-T"
# if strand == -, complement both alleles
# ----------------------------
nba_to_cpra <- function(coord, base_change) {
  if (is.na(coord) || coord == "" || is.na(base_change) || base_change == "") {
    return(NA_character_)
  }
  
  chr <- str_match(coord, "^chr([0-9XYM]+):")[,2]
  pos <- str_match(coord, "^chr[0-9XYM]+:([0-9]+)")[,2]
  strand <- str_match(coord, "\\(([+-])\\)")[,2]
  
  bc <- str_split(base_change, "-", simplify = TRUE)
  if (ncol(bc) < 2) return(NA_character_)
  
  ref <- toupper(trimws(bc[1]))
  alt <- toupper(trimws(bc[2]))
  
  if (is.na(chr) || is.na(pos) || is.na(strand)) return(NA_character_)
  
  if (strand == "-") {
    ref <- comp_base(ref)
    alt <- comp_base(alt)
  }
  
  paste(chr, pos, ref, alt, sep = ":")
}

# ----------------------------
# 4) Genoox string -> CPRA
# example:
# chr1-7969392 GAAT>G  ->  1:7969392:GAAT:G
# ----------------------------
genoox_to_cpra <- function(x) {
  if (is.na(x) || x == "") return(NA_character_)
  
  m <- str_match(x, "^chr([0-9XYM]+)-([0-9]+)\\s+([^>]+)>(.+)$")
  chr <- m[,2]
  pos <- m[,3]
  ref <- m[,4]
  alt <- m[,5]
  
  if (is.na(chr) || is.na(pos) || is.na(ref) || is.na(alt)) return(NA_character_)
  
  paste(chr, pos, ref, alt, sep = ":")
}

# ----------------------------
# 5) build final comparison df
# ----------------------------
cpra_df <- toannot %>%
  transmute(
    unique_variants,
    CPRA_CLINVAR = sapply(.data$`Canonical SPDI`, spdi_to_cpra),
    CPRA_NBA     = mapply(
      nba_to_cpra,
      .data$`Genomic.coordinates..GRCh38.`,
      .data$`Base.change`,
      USE.NAMES = FALSE
    ),
    CPRA_GENOOX  = sapply(.data$genoox, genoox_to_cpra)
  ) %>%
  distinct(unique_variants, .keep_all = TRUE)



cpra_final <- cpra_df %>%
  mutate(
    CPRA = coalesce(CPRA_GENOOX, CPRA_CLINVAR, CPRA_NBA)
  ) %>%
  select(unique_variants, CPRA, CPRA_GENOOX, CPRA_CLINVAR, CPRA_NBA)



missing <- cpra_final %>% filter(is.na(CPRA))
missing


# ==========================================================
# 5. Annovar Preparation
# ==========================================================


##################################################### ANNOVAR PREP

annovar_input <- cpra_final %>%
  filter(!is.na(CPRA)) %>%
  separate(CPRA, into = c("chr", "pos", "ref", "alt"), sep = ":", remove = FALSE) %>%
  mutate(
    pos = as.numeric(pos),
    start = pos,
    end = pos + nchar(ref) - 1
  ) %>%
  select(chr, start, end, ref, alt, unique_variants, CPRA)



vcf_df <- cpra_final %>%
  filter(!is.na(CPRA)) %>%
  separate(CPRA, into = c("CHROM", "POS", "REF", "ALT"), sep = ":", remove = FALSE) %>%
  mutate(
    CHROM = paste0("chr", CHROM),
    POS = as.integer(POS),
    ID = unique_variants,
    REF = ifelse(is.na(REF) | REF == "", "-", REF),
    ALT = ifelse(is.na(ALT) | ALT == "", "-", ALT),
    QUAL = ".",
    FILTER = ".",
    INFO = "."
  ) %>%
  select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO)

contigs <- sort(unique(vcf_df$CHROM))

vcf_path <- "variants_raw_DIGENIC.vcf"

con <- file(vcf_path, "w")
writeLines("##fileformat=VCFv4.2", con)
writeLines("##source=cpra_from_manual_merge", con)
writeLines(paste0("##contig=<ID=", contigs, ">"), con)
writeLines("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", con)
close(con)

write.table(
  vcf_df,
  file = vcf_path,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  append = TRUE
)

############## Next section of code is executable in Bash


#bcftools norm -f hg38.fa -Oz -o variants_norm_DIGENIC.vcf.gz variants_raw_DIGENIC.vcf

#perl convert2annovar.pl \
#-format vcf4 \
#-allsample \
#-withfreq \
#variants_norm_DIGENIC.vcf.gz \
#> variants_norm_DIGENIC.avinput



#perl table_annovar.pl \
#variants_norm_DIGENIC.avinput \
#/utils/annovar/humandb/ \
#-buildver hg38 \
#-out variants_annotated_DIGENIC \
#-remove \
#-protocol refGene,clinvar_20250721,dbnsfp47a \
#-operation g,f,f \
#-nastring . \
#-polish \
#-xreffile /utils/annovar/example/gene_xref.txt


#bcftools view variants_norm_DIGENIC.vcf.gz | grep -v "##" > variants_DIGENIC_noralised_vcf_IDs.tsv



# ==========================================================
# 6. Variant Filtering and Cleaning
# ==========================================================


################################################## Adding annot info
norm_Ids <- read.csv("variants_DIGENIC_noralised_vcf_IDs.tsv", sep = "\t", header = TRUE)


norm_Ids$unique_variants <- norm_Ids$ID
norm_Ids$CPRA <- paste0(norm_Ids$X.CHROM,":", norm_Ids$POS, ":", norm_Ids$REF, ":", norm_Ids$ALT)


annovar_raw <- read.csv("variants_annotated_DIGENIC.hg38_multianno.txt", sep = "\t", header = TRUE)

annovar <- bind_cols(annovar_raw,norm_Ids)
colnames(annovar)

cols_to_keep <- c(
  "Chr", "Start", "End",
  "Ref", "Alt",
  "Func.refGene", "Gene.refGene", "GeneDetail.refGene",
  "ExonicFunc.refGene", "AAChange.refGene",
  "unique_variants", "CPRA","CLNALLELEID",	"CLNDN",	"CLNDISDB",	"CLNREVSTAT",	"CLNSIG"
)

annovar <- annovar[,cols_to_keep]

final_annot <- annovar

write.csv(final_annot, 
          "FINAL_ANNOTATED_VARIANTS_DIGENIC.csv", 
          row.names = FALSE)

colnames(final_annot)



table(final_annot$Func.refGene)

table(final_annot$ExonicFunc.refGene)


########################################################## Harmonising variants across CH Px table

####### checking which variants have more that one refseq ID

dup_vars <- final_annot[duplicated(final_annot$CPRA) | duplicated(final_annot$CPRA, fromLast = TRUE), 
                        c("unique_variants", "CPRA", "Gene.refGene", "AAChange.refGene")]

dup_vars <- dup_vars[order(dup_vars$CPRA), ]
dup_vars
cpra_counts <- as.data.frame(table(final_annot$CPRA), stringsAsFactors = FALSE)
cpra_counts <- cpra_counts[cpra_counts$Freq > 1, ]
cpra_counts

######### no variants with more than one refseq ID

########################## sanity check
## FILTERED OUT BENING VARIANTS

benign_vars <- final_annot %>%
  filter(str_detect(CLNSIG, regex("benign", ignore_case = TRUE))) %>%
  pull(unique_variants) %>%
  unique()

di_filtered <- di_raw %>%
  filter(
    !var_1_g %in% benign_vars,
    !var_2_g %in% benign_vars
  )

##Checking homozygotes
di_filtered <- di_filtered %>%
  filter(var_1_g != var_2_g)  # remove homozygotes

## Excluding cases without enough genetic data ( no annotation available)
fix_sv <- function(x) {
  x %>%
    str_replace_all("–", "-") %>% 
    str_replace_all("Ex(\\d+)_([0-9P]+)", "Ex\\1-\\2")
}

is_rearrangement <- function(x) {
  grepl("_(Ex|ex).*(del|dup|inv|partial)", x)
}

di_filtered <- di_filtered %>%
  mutate(
    var_1_g = fix_sv(var_1_g),
    var_2_g = fix_sv(var_2_g)
  )

annot_vars <- unique(final_annot$unique_variants)

di_filtered <- di_filtered %>%
  filter(
    (var_1_g %in% annot_vars | is_rearrangement(var_1_g)),
    (var_2_g %in% annot_vars | is_rearrangement(var_2_g))
  )

################### getting rid of cases where one of the variants is associated eith Autosomal dominant PD

####
final_annot %>%
  filter(grepl("dominant", CLNDN, ignore.case = TRUE)) %>%
  pull(Gene.refGene)

table(final_annot$Gene.refGene)

di_filtered_dom <- di_filtered %>%
  filter( Gene_1 != "LRRK2" & Gene_2 != "LRRK2")

###########################


vars_all <- unique(c(as.character(di_filtered$var_1_g), as.character(di_filtered$var_2_g)))
vars_all <- vars_all[!is.na(vars_all)]

is_rearrangement <- grepl("_(Ex|ex).*(del|dup|inv|partial)", vars_all)

vars_rearr <- sort(vars_all[is_rearrangement])
vars_seq   <- sort(vars_all[!is_rearrangement])



vars_all_clean <-  tibble(full = vars_all) %>%
  separate(full, into = c("gene", "variant"), sep = "_", extra = "merge")

SVs_clean <-  tibble(full = vars_rearr) %>%
  separate(full, into = c("gene", "variant"), sep = "_", extra = "merge")

SNVs_clean <-  tibble(full = vars_seq) %>%
  separate(full, into = c("gene", "variant"), sep = "_", extra = "merge")



SNVs_clean$unique_variants <- paste0(SNVs_clean$gene,"_", SNVs_clean$variant)
SNVs_clean <- merge(SNVs_clean, final_annot, by = "unique_variants")
SVs_clean$unique_variants <- paste0(SVs_clean$gene,"_", SVs_clean$variant)



missing_seq <- setdiff(as.character(vars_seq), as.character(final_annot$unique_variants))
missing_seq



write.csv(di_filtered, 
          "DI_CLEANLIST.csv", 
          row.names = FALSE)
write.csv(SNVs_clean, 
          "DI_SNVs.csv", 
          row.names = FALSE)

write.csv(SVs_clean, 
          "DI_SVs.csv", 
          row.names = FALSE)


# ==========================================================
# 7. Manuscript Statistics
# ==========================================================

length(vars_all_clean$variant)


## exonic function classification per SNV and small indel (no SVs)
length(SNVs_clean$unique_variants)
table(SNVs_clean$ExonicFunc.refGene)

### Number of rearrangements
length(SVs_clean$variant)


table(str_extract(unique(SVs_clean$variant), "del|dup|inv")) ##types or rearrs in PRKN



############################################################# LRRK2-prkn
lrkk2_prkn <- di_filtered %>%
  mutate(
    gene_A = pmin(Gene_1, Gene_2),
    gene_B = pmax(Gene_1, Gene_2)
  ) %>%
  filter(gene_A == "LRRK2", gene_B == "PRKN")

lrkk2_prkn_varcount <- lrkk2_prkn %>%
  select(var_1_g, var_2_g) %>%
  pivot_longer(cols = everything(), values_to = "variant") %>%
  count(variant, sort = TRUE)

lrkk2_prkn_varcount

repeated_combos <- lrkk2_prkn %>%
  mutate(
    var_A = pmin(var_1_g, var_2_g),
    var_B = pmax(var_1_g, var_2_g)
  ) %>%
  count(var_A, var_B, sort = TRUE) %>%
  filter(n > 1)

repeated_combos
############################################### check pairs for any given variant
get_variant_combos <- function(df, variant) {
  df %>%
    filter(var_1_g == variant | var_2_g == variant) %>%
    mutate(
      partner = ifelse(var_1_g == variant, var_2_g, var_1_g)
    ) %>%
    count(partner, sort = TRUE)
}

get_variant_combos(lrkk2_prkn, "PRKN_p.Met192Val")
lrkk2_prkn <- di_filtered %>%
  mutate(
    gene_A = pmin(Gene_1, Gene_2),
    gene_B = pmax(Gene_1, Gene_2)
  ) %>%
  filter(gene_A == "LRRK2", gene_B == "PRKN")

lrkk2_prkn_varcount <- lrkk2_prkn %>%
  select(var_1_g, var_2_g) %>%
  pivot_longer(cols = everything(), values_to = "variant") %>%
  count(variant, sort = TRUE)

lrkk2_prkn_varcount

repeated_combos <- lrkk2_prkn %>%
  mutate(
    var_A = pmin(var_1_g, var_2_g),
    var_B = pmax(var_1_g, var_2_g)
  ) %>%
  count(var_A, var_B, sort = TRUE) %>%
  filter(n > 1)

repeated_combos
############################################### check pairs for any given variant
get_variant_combos <- function(df, variant) {
  df %>%
    filter(var_1_g == variant | var_2_g == variant) %>%
    mutate(
      partner = ifelse(var_1_g == variant, var_2_g, var_1_g)
    ) %>%
    count(partner, sort = TRUE)
}

get_variant_combos(lrkk2_prkn, "PRKN_p.Met192Val")



#############################################population stratification
lrkk2_prkn <- lrkk2_prkn %>%
  mutate(
    Population_clean = str_to_lower(Population),
    Population_continent = case_when(
      
      # EUROPE
      str_detect(Population_clean, "europe|caucasian|french|german|irish|italian|italy|polish|serbian|swiss|uk|ukrainian|czech|netherlands|polish|finnland|france|belgian|spain|norway|luxembourg") ~ "EUR",
      
      #Middle East
      str_detect(Population_clean, "cyprus|turkish|yemen") ~ "MDE",
      
      # EAST ASIA
      str_detect(Population_clean, "china|chinese|japan|japanese|korea|korean|taiwan|vietnam|filipino") ~ "EAS",
      
      # AFRICA
      str_detect(Population_clean, "african|tunisia") ~ "AFR",
      
      # LATIN AMERICA / ADMIXED AMERICAN
      str_detect(Population_clean, "mexico|mexican|brazil|peru") ~ "AMR",
      
      # NORTH AMERICA (ambiguous, likely EUR )
      str_detect(Population_clean, "north american") ~ "EUR",
      
      TRUE ~ "Other"
    )
  )



table(lrkk2_prkn$Population_continent)







############################################################# pink1-prkn
pink1_prkn <- di_filtered %>%
  mutate(
    gene_A = pmin(Gene_1, Gene_2),
    gene_B = pmax(Gene_1, Gene_2)
  ) %>%
  filter(gene_A == "PINK1", gene_B == "PRKN")

pink1_prkn_varcount <- pink1_prkn %>%
  select(var_1_g, var_2_g) %>%
  pivot_longer(cols = everything(), values_to = "variant") %>%
  count(variant, sort = TRUE)

pink1_prkn_varcount

repeated_combos <- pink1_prkn %>%
  mutate(
    var_A = pmin(var_1_g, var_2_g),
    var_B = pmax(var_1_g, var_2_g)
  ) %>%
  count(var_A, var_B, sort = TRUE) %>%
  filter(n > 1)

repeated_combos

############################################### check pairs for any given variant
get_variant_combos <- function(df, variant) {
  df %>%
    filter(var_1_g == variant | var_2_g == variant) %>%
    mutate(
      partner = ifelse(var_1_g == variant, var_2_g, var_1_g)
    ) %>%
    count(partner, sort = TRUE)
}

get_variant_combos(pink1_prkn, "PRKN_p.Met192Val")

#############################################population stratification
pink1_prkn <- pink1_prkn %>%
  mutate(
    Population_clean = str_to_lower(Population),
    Population_continent = case_when(
      
      # EUROPE
      str_detect(Population_clean, "europe|caucasian|french|german|irish|italian|italy|polish|serbian|swiss|uk|ukrainian|czech|netherlands|polish|finnland|france|belgian|spain|norway|luxembourg") ~ "EUR",
      
      #Middle East
      str_detect(Population_clean, "cyprus|turkish|yemen") ~ "MDE",
      
      # EAST ASIA
      str_detect(Population_clean, "china|chinese|japan|japanese|korea|korean|taiwan|vietnam|filipino") ~ "EAS",
      
      # AFRICA
      str_detect(Population_clean, "african|tunisia") ~ "AFR",
      
      # LATIN AMERICA / ADMIXED AMERICAN
      str_detect(Population_clean, "mexico|mexican|brazil|peru") ~ "AMR",
      
      # NORTH AMERICA (ambiguous, likely EUR )
      str_detect(Population_clean, "north american") ~ "EUR",
      
      TRUE ~ "Other"
    )
  )



table(pink1_prkn$Population_continent)




############################################################# pink1-dj-1
pink1_park7 <- di_filtered %>%
  mutate(
    gene_A = pmin(Gene_1, Gene_2),
    gene_B = pmax(Gene_1, Gene_2)
  ) %>%
  filter(gene_A == "PARK7", gene_B == "PINK1")

pink1_park7_varcount <- pink1_park7 %>%
  select(var_1_g, var_2_g) %>%
  pivot_longer(cols = everything(), values_to = "variant") %>%
  count(variant, sort = TRUE)

pink1_park7_varcount

repeated_combos <- pink1_park7 %>%
  mutate(
    var_A = pmin(var_1_g, var_2_g),
    var_B = pmax(var_1_g, var_2_g)
  ) %>%
  count(var_A, var_B, sort = TRUE) %>%
  filter(n > 1)

repeated_combos



############################################################# OTHERS
others <- di_filtered %>%
  mutate(
    gene_A = pmin(Gene_1, Gene_2),
    gene_B = pmax(Gene_1, Gene_2)
  ) %>%
  filter(
    !(gene_A == "LRRK2" & gene_B == "PRKN"),
    !(gene_A == "PINK1" & gene_B == "PRKN"),
    !(gene_A == "PARK7" & gene_B == "PINK1")
  )

others_varcount <- others %>%
  select(var_1_g, var_2_g) %>%
  pivot_longer(cols = everything(), values_to = "variant") %>%
  count(variant, sort = TRUE)

others_varcount

repeated_combos <- others %>%
  mutate(
    var_A = pmin(var_1_g, var_2_g),
    var_B = pmax(var_1_g, var_2_g)
  ) %>%
  count(var_A, var_B, sort = TRUE) %>%
  filter(n > 1)

repeated_combos

############################################### check pairs for any given variant
get_variant_combos <- function(df, variant) {
  df %>%
    filter(var_1_g == variant | var_2_g == variant) %>%
    mutate(
      partner = ifelse(var_1_g == variant, var_2_g, var_1_g)
    ) %>%
    count(partner, sort = TRUE)
}

get_variant_combos(others, "PRKN_p.Met192Val")

#############################################population stratification
others <- others %>%
  mutate(
    Population_clean = str_to_lower(Population),
    Population_continent = case_when(
      
      # EUROPE
      str_detect(Population_clean, "europe|caucasian|french|german|irish|italian|italy|polish|serbian|swiss|uk|ukrainian|czech|netherlands|polish|finnland|france|belgian|spain|norway|luxembourg") ~ "EUR",
      
      #Middle East
      str_detect(Population_clean, "cyprus|turkish|yemen") ~ "MDE",
      
      # EAST ASIA
      str_detect(Population_clean, "china|chinese|japan|japanese|korea|korean|taiwan|vietnam|filipino") ~ "EAS",
      
      # AFRICA
      str_detect(Population_clean, "african|tunisia") ~ "AFR",
      
      # LATIN AMERICA / ADMIXED AMERICAN
      str_detect(Population_clean, "mexico|mexican|brazil|peru") ~ "AMR",
      
      # NORTH AMERICA (ambiguous, likely EUR )
      str_detect(Population_clean, "north american") ~ "EUR",
      
      TRUE ~ "Other"
    )
  )



table(others$Population_continent)


dosi_digenic <- as.data.frame(unique(di_filtered$DOI))
colnames(dosi_digenic)[colnames(dosi_digenic) == "unique(di_filtered$DOI)"] <- "DOI"

overlap_dois <- dplyr::inner_join(inccluded_papers, dosi_digenic, by = "DOI")

query_digenic <- read.csv("csv-digenicAND-set_2026.csv", header = TRUE)

included_di <- merge(query_digenic, dosi_digenic, by = "DOI")

write.csv(included_di, "included_papers_DIGENIC_n9.csv", row.names = FALSE)



##### redo DI cleanlist

rating_di <- read.csv("included_papers_DIGENIC_RATING_n9.csv", header = TRUE)

rating_di <- rating_di[,c("DOI", "Rating")]
table(rating_di$Rating)

di_filtered_rating <- merge(rating_di, di_filtered, by = "DOI", all.x = TRUE)

di_filtered_rating <- di_filtered_rating[order(di_filtered_rating$IID), ]

write.csv(di_filtered_rating, 
          "DI_CLEANLIST_RATING_APR24.csv", 
          row.names = FALSE)
