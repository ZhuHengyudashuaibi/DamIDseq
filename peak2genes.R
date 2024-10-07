# Check pakage requirement
check_pkg <- function(pkg){
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE, quietly = T)
  } else {
    library(pkg, character.only = TRUE, quietly = T)
  }
}

check_pkg_bioconductor <- function(pkg){
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE, quietly = T)
  } else {
    library(pkg, character.only = TRUE, quietly = T)
  }
  
}


# Define function
peak_overlap <- function(sig_peak, combine_mode = "max"){
  all_peak <- read.table(sig_peak, header = F, sep = "\t")
  i <- 0
  for (chr in unique(all_peak[[1]])) {
    each_chr <- all_peak[all_peak[[1]] == chr,]
    for (pk in 1:nrow(each_chr)) {
      each_pk <- each_chr[pk, ]
      overlap <- each_chr %>% filter(V2 >= each_pk$V2, V3 <= each_pk$V3) %>% as.data.frame()
      if (combine_mode == "max") {
        start <- min(overlap$V2)
        end <- max(overlap$V3)
      }else {
        start <- max(overlap$V2)
        end <- min(overlap$V3)
      }
      score <- max(overlap$V4)
      fdr <- min(overlap$V8)
      if (i == 0) {
        unique_peak <- data.frame(chr = chr, start = start, end = end, score = score, fdr = fdr)
      }else {
        unique_peak <- rbind(unique_peak, 
                             data.frame(chr = chr, start = start, end = end, score = score, fdr = fdr))
      }
      i <- i + 1
    }
  }
  unique_peak <- unique_peak %>% arrange(chr, start, end)
  return(unique_peak)
}

peak2gene <- function(peak, gxf, gene_pad = 1000, up_pad = 6000, down_pad = 1000, mode = "promoter", anno_type = "gene"){
  gxf_data <- import(gxf) %>% as.data.frame()
  gene_info <- gxf_data %>% filter(type == anno_type)
  if (nrow(gene_info) == 0) {
    stop("Error: wrong annotation type or wrong gxf file, please check!")
  }
  i <- 0
  for (chr in unique(peak$chr)) {
    each_chrpk <- peak[peak$chr == chr,]
    for (pk in 1:nrow(each_chrpk)) {
      if (mode == "promoter") {
        pos_gene <- gene_info %>% filter(seqnames == chr, strand == "+") %>%
                                  filter(start <= each_chrpk[pk,]$end + gene_pad, 
                                         start >= each_chrpk[pk,]$start - gene_pad) %>% 
                                  select(seqnames, Name)
        neg_gene <- gene_info %>% filter(seqnames == chr, strand == "-") %>%
                                  filter(end <= each_chrpk[pk,]$end + gene_pad, 
                                          end >= each_chrpk[pk,]$start - gene_pad)%>% 
                                  select(seqnames, Name)
        enriched <- rbind(pos_gene, neg_gene)
      }else if (mode == "loose"){
        enriched <- gene_info %>% filter(seqnames == chr) %>%
                                  filter(start <= each_chrpk[pk,]$end + gene_pad, 
                                         end >= each_chrpk[pk,]$start - gene_pad) %>% 
                                  select(seqnames, Name)
      }else if (mode == "center") {
        # Center mode reference: Katsanos and Barkoulas, Sci. Adv. 8, eabk3141 (2022)
        center_loc <- (each_chrpk[pk,]$end + each_chrpk[pk,]$start)/2
        
        pos_gene <- gene_info %>% filter(seqnames == chr, strand == "+") %>%
                                  filter(start <= center_loc + up_pad, 
                                         end >= center_loc - down_pad) %>% 
                                  select(seqnames, Name)
        neg_gene <- gene_info %>% filter(seqnames == chr, strand == "-") %>%
                                  filter(start <= center_loc + down_pad, 
                                         end >= center_loc - up_pad)%>% 
                                  select(seqnames, Name)
        enriched <- rbind(pos_gene, neg_gene)
      }else {
        stop("Unsupported assignment mode. Please check.")
      }
      
      if (nrow(enriched) != 0) {
        enriched$score <- each_chrpk[pk,]$score
        enriched$fdr <- each_chrpk[pk,]$fdr
        
        # Assign genes
        types <- c()
        for (each_gene in enriched$Name) {
          strand <- filter(gene_info, Name == each_gene)$strand
          
          if (strand == "-") {
            gene.st <- filter(gene_info, Name == each_gene)$end
            gene.ed <- filter(gene_info, Name == each_gene)$start
            pk.st <- each_chrpk[pk,]$start
            pk.ed <- each_chrpk[pk,]$end
            if (pk.st >= gene.st) {
              type <- "Upstream"
            }else if (pk.st < gene.st & pk.ed > gene.st) {
              type <- "Overlap with TSS"
            }else if (pk.ed <= gene.st & pk.st >= gene.ed) {
              type <- "Peak in gene"
              
              # is.exon
              gene_ID <- gxf_data %>% filter(Name == each_gene | grepl("Y47G6A.28", Name))
              gene_ID <- gxf_data %>% filter(Name == each_gene)
              gene_ID <- gxf_data %>% filter(grepl("Y47G6A.28", transcript_id))
              gene_ID <- gsub("^[^:]+:", "", gene_ID)
              
              
            }else if (pk.st < gene.ed & pk.ed > gene.ed) {
              type <- "Overlap with TES"
            }else if (pk.ed <= gene.ed) {
              type <- "Downstream"
            }else if (pk.st <= gene.ed & pk.ed >= gene.st) {
              type <- "Gene in peak"
            }else {
              stop("Please check you code.")
            }
          }else if (strand == "+") {
            gene.st <- filter(gene_info, Name == each_gene)$start
            gene.ed <- filter(gene_info, Name == each_gene)$end
            pk.st <- each_chrpk[pk,]$start
            pk.ed <- each_chrpk[pk,]$end
            if (pk.ed <= gene.st) {
              type <- "Upstream"
            }else if (pk.st < gene.st & pk.ed > gene.st) {
              type <- "Overlap with TSS"
            }else if (pk.ed <= gene.ed & pk.st >= gene.st) {
              type <- "Peak in gene"
            }else if (pk.st < gene.ed & pk.ed > gene.ed) {
              type <- "Overlap with TES"
            }else if (pk.st >= gene.ed) {
              type <- "Downstream"
            }else if (pk.st <= gene.st & pk.ed >= gene.ed) {
              type <- "Gene in peak"
            }else {
              stop("Please check you code.")
            }
          }
          types <- c(types, type)
        }
        enriched$type <- types
        
        if (i == 0) {
          enriched_all <- enriched
        }else {
          enriched_all <- rbind(enriched_all, enriched)
        }
        i <- i + 1
        cat(sprintf("\rNow processing: %d", i))
      }
    }
  }
  
  enriched_all <- enriched_all %>% group_by(Name) %>% 
                  arrange(desc(score), fdr) %>% 
                  filter(row_number() == 1) %>% ungroup()
  
  
  return(enriched_all)
}

# resolve arguments
args <- commandArgs(trailingOnly = TRUE)

show_helpdoc <- function(){
  cat("\nusage: Rscript peak2gene.R 
      --gxf=<Path to annotation file> \\\
      --significant_peaks=<Path to significant peaks> \\\
      --combine_mode=<max/随便什么> \\\
      --anno_type=gene \\\
      --assign_mode=<promoter/loose/center> \\\
      --gene_pad=1000 \\\
      --gene_up_pad=6000 \\\
      --gene_down_pad=1000\n")
  
  cat("\n------------Arguments-------------\n")
  cat("\n--gxf: gxf file path, should be gff3 or gtf. (Default: ~/zhuhengyu/DCL1125_1111_damID/damID/damid_data/peak/peak_analysis.avg.2024-08-20.11-20-35/avg-FDR0.05allpeaks.gff)\n")
  cat("\n--significant_peaks: path to significant peaks data, recommended to be generated by find_peaks.py.\n")
  cat("\n--combine_mode: mode of overlapped peaks combination. If this argument is max, the most upstream coordination and the most downstream coordination of multiple overlapped peaks will be considered to be the start/end site of the combined peak. Otherwise, only the shared region of overlapped peaks will be used.(Default is max)\n")
  cat("\n--assign_mode: mode of peak assignment, should be any one of promoter, loose or center (default is center). \n
       \npromoter: assign a peak to a gene when any part of the peak overlap with the transcription start site ± gene_pad. \n
        \nloose: assign a peak to a gene when any part of the peak overlap with the region between transcription start site - gene_pad and transcription end site + gene_pad. \n
        \ncenter: assign a peak to a gene only when the center coordination is within the region between TSS - gene_up_pad and TES + gene_down_pad. \n")
  cat("\n--gene_pad: only useful when assign_mode=promoter/loose. Length of chromosome that should be considered as the functional regulatory region around the TSS/TES (default is 100). \n")
  cat("\n--gene_up_pad: only useful when assign_mode=center. Length of chromosome that should be considered as the functional regulatroy region upstream of TSS (default is 6000). \n")
  cat("\n--gene_down_pad: only useful when assign_mode=center. Length of chromosome that should be considered as the functional regulatroy region downstream of TES (default is 1000). \n")
  cat("\n--anno_type: a character that represent the whole gene within the column 7 of the gxf file (default is gene). \n")
  cat("\n--help: show help document.\n")
}

resolve_arg <- function(argument_type, args=commandArgs(trailingOnly = TRUE), default){
  target <- args[grepl(argument_type, args)]
  if (length(target) == 0) {
    value <- default
  }else {
    value <- gsub(paste0("^", argument_type, "="), "", target)
  }
  return(value)
}

if (length(args) == 0 | length(args[grepl("--help", args)]) == 1) {
  show_helpdoc()
}else if (length(args) == 1 & length(args[grepl("--significant_peaks", args)]) == 1) {
  gxf <- "~/reference_genome/gxf/Caenorhabditis_elegans.WBcel235.112.gff3"
  combine_mode <- "max"
  assign_mode <- "center"
  gene_pad <- 1000
  gene_up_pad <- 6000
  gene_down_pad <- 1000
  anno_type <- "gene"
  
  check_pkg_bioconductor("rtracklayer")
  check_pkg("dplyr")
  
  significant_peaks <- resolve_arg("--significant_peaks", args)
  peak <- peak_overlap(sig_peak = significant_peaks)
  assignment <- peak2gene(peak = peak, gxf = gxf, gene_pad = gene_pad, up_pad = gene_up_pad, down_pad = gene_down_pad, mode = assign_mode, anno_type = anno_type)
  write.csv(assignment, row.names = F, file = paste0("./peak_assignment_", anno_type, ".csv"))
  
}else if (length(args[grepl("--significant_peaks", args)]) == 0) {
 stop("Error: need to parse the necessary argument '--significant_peaks'") 
}else {
  
  check_pkg_bioconductor("rtracklayer")
  check_pkg("dplyr")
  
  gxf <- resolve_arg("--gxf", args, "~/reference_genome/gxf/Caenorhabditis_elegans.WBcel235.112.gff3")
  significant_peaks <- resolve_arg("--significant_peaks", args)
  combine_mode <- resolve_arg("--combine_mode", args, "center")
  anno_type <- resolve_arg("--anno_type", args, "gene")
  assign_mode <- resolve_arg("--assign_mode", args, "max")
  gene_pad <- as.numeric(resolve_arg("--gene_pad", args, "1000"))
  gene_up_pad <- as.numeric(resolve_arg("--gene_up_pad", args, "6000"))
  gene_down_pad <- as.numeric(resolve_arg("--gene_down_pad", args, "1000"))
  
  # Peak assignment
  
  peak <- peak_overlap(sig_peak = significant_peaks)
  assignment <- peak2gene(peak = peak, gxf = gxf, gene_pad = gene_pad, up_pad = gene_up_pad, down_pad = gene_down_pad, mode = assign_mode, anno_type = anno_type)
  write.csv(assignment, row.names = F, file = paste0("./peak_assignment_", anno_type, ".csv"))
  
}