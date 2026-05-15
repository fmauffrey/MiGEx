reads_quality_plot <- function(fastp_path){
  # Reads fastp log and return histogram plot on the number of filtered reads
  log_text <- readChar(fastp_path, file.info(fastp_path)$size)
  
  ##### First plot - number of filtered reads
  text_reads_number <- str_extract_all(log_text, "Read(.|\n)*?\\d\n")
  
  reads_data <- data.frame(
    read_file = str_extract(unlist(text_reads_number), "^Read\\d"),
    state = str_extract(unlist(text_reads_number), "before|after"),
    value = as.numeric(str_extract(unlist(text_reads_number), "\\d+\n")))
  
  reads_data <- reads_data %>%
    mutate(read_file = case_when(
      read_file == "Read1" ~ "R1",
      read_file == "Read2" ~ "R2"
    )) %>%
    mutate(state = factor(state, levels=c("before", "after")))
  
  reads_plot <- ggplot(reads_data, aes(x=read_file, y=value, fill=state))+
    geom_histogram(stat = "identity", position = 'dodge')+
    theme_classic()+
    labs(x="", fill="", y="Reads number")+
    scale_fill_manual(values = c("before" = "purple3",
                                 "after" = "orange3"))+
    geom_text(
      aes(label = value),  
      position = position_dodge(width = 0.9), 
      vjust = 1.5,
      color = "white",
      size = 3.5)
  
  #### Second plot - quality
  Q20_text <- str_extract_all(log_text, "Read\\d after(.|\n)*?Q20.*?%")
  Q30_text <- str_extract_all(log_text, "Read\\d after(.|\n)*?Q30.*?%")
  Q40_text <- str_extract_all(log_text, "Read\\d after(.|\n)*?Q40.*?%")
  
  Q20_data <- data.frame(
    read_file = str_extract(unlist(Q20_text), "^Read\\d"),
    qual = rep("Q20", 1),
    value = as.numeric(str_extract(unlist(Q20_text), "\\d+\\.\\d+(?=%$)")))
  
  Q30_data <- data.frame(
    read_file = str_extract(unlist(Q30_text), "^Read\\d"),
    qual = rep("Q30", 1),
    value = as.numeric(str_extract(unlist(Q30_text), "\\d+\\.\\d+(?=%$)")))
  
  Q40_data <- data.frame(
    read_file = str_extract(unlist(Q40_text), "^Read\\d"),
    qual = rep("Q40", 1),
    value = as.numeric(str_extract(unlist(Q40_text), "\\d+\\.\\d+(?=%$)")))
  
  quality_data <- rbind.data.frame(Q20_data, Q30_data, Q40_data) %>%
    mutate(read_file = case_when(
      read_file == "Read1" ~ "R1",
      read_file == "Read2" ~ "R2"
    )) %>%
    mutate(value = round(value, 1))
  
  quality_plot_R1 <- ggplot(quality_data[quality_data$read_file=="R1",],
                            aes(x=value, y=qual, fill=qual))+
    geom_histogram(stat = "identity")+
    theme_classic()+
    labs(x="Reads number (%)", y="Quality")+
    scale_fill_manual(values = c("Q20" = "#2F6280",
                                 "Q30" = "#5EA1B0",
                                 "Q40" = "#B0DFDD"))+
    geom_text(label=quality_data[quality_data$read_file=="R1",3],
              hjust = 1,
              size = 5,
              color = "white")+
    theme(legend.position="none")
  
  quality_plot_R2 <- ggplot(quality_data[quality_data$read_file=="R2",],
                            aes(x=value, y=qual, fill=qual))+
    geom_histogram(stat = "identity")+
    theme_classic()+
    labs(x="Reads number (%)", y="Quality")+
    scale_fill_manual(values = c("Q20" = "#2F6280",
                                 "Q30" = "#5EA1B0",
                                 "Q40" = "#B0DFDD"))+
    geom_text(label=quality_data[quality_data$read_file=="R1",3],
              hjust = 1,
              size = 5,
              color = "white")+
    theme(legend.position="none")
  
  second_panel <- ggarrange(quality_plot_R1, quality_plot_R2, ncol=1, labels = c("R1", "R2"))
  
  ggarrange(reads_plot, second_panel)
}

contigs_assembly_plot <- function(quast_path, coverage_path, pipeline_path){
  # Plot for assembly quality metrics
  
  # Load files
  quast <- read.csv(quast_path, sep="\t")
  coverage <- read.csv(coverage_path, sep="\t")
  
  #### Plot contigs number (left panel)
  contigs_data <- quast[1:6,] %>%
    mutate(Assembly = factor(c(">= 0kb", ">= 1kb", ">= 5kb", ">= 10kb", ">= 25kb", ">= 50kb"), 
                             levels=c(">= 0kb", ">= 1kb", ">= 5kb", ">= 10kb", ">= 25kb", ">= 50kb")))
  
  contigs_number_plot <- ggplot(contigs_data, aes(x=Assembly, y=assembly))+
    geom_histogram(stat = "identity")+
    theme_minimal()+
    labs(x="Contig size", y="Contigs number")+
    geom_text(label=contigs_data$assembly,
              vjust=1.2,
              color="white",
              size=6)
  
  #### Generate metrics plot (right panel)
  # Get metrics
  total_length <- as.numeric(quast[quast$Assembly=="Total length",2])
  GC <- as.numeric(quast[quast$Assembly=="GC (%)",2])
  mean_depth <- round(weighted.mean(x=coverage$meandepth, coverage$endpos),2)
  
  total_length_text <- paste0(comma(total_length), "bp")
  GC_text <- paste0(GC, "%")
  mean_depth_text <- paste0(mean_depth, "x")
  
  # Load icons and convert them 
  img_total_length <- readPNG(paste0(pipeline_path, "/img/icon_total_length.png"))
  img_GC <- readPNG(paste0(pipeline_path, "/img/icon_GC.png"))
  img_mean_depth <- readPNG(paste0(pipeline_path, "/img/icon_mean_depth.png"))
  
  grob_total_length <- rasterGrob(img_total_length, width = unit(70, "points"),
                                  height = unit(50, "points"),
                                  interpolate = TRUE)
  grob_GC <- rasterGrob(img_GC, width = unit(50, "points"),
                        height = unit(50, "points"),
                        interpolate = TRUE)
  grob_mean_depth <- rasterGrob(img_mean_depth, width = unit(70, "points"),
                                height = unit(50, "points"),
                                interpolate = TRUE)

  # Icons column
  icons_col <- plot_grid(
    grob_total_length,
    grob_GC,
    grob_mean_depth,
    ncol = 1
  )
  
  make_text <- function(text){
    ggdraw()+
      draw_label(text, x = 0, hjust = 0, size = 20)+
      theme(plot.margin = margin(0, 0, 0, 0))
  }
  
  text_col <- plot_grid(
    make_text(total_length_text),
    make_text(GC_text),
    make_text(mean_depth_text),
    ncol = 1
  )
  
  left_panel <- plot_grid(
    icons_col,
    text_col,
    ncol = 2,
    rel_widths = c(1,1)
  )
  
  final_plot <- plot_grid(
    left_panel,
    contigs_number_plot,
    ncol = 2,
    rel_widths = c(1, 2)
  )
  
  final_plot
}

identification_table <- function(mash_path, mash_details_path){
  # Extract relevant information from the JSON file and return a table
  
  # Top hits identification
  mash_details_data <- fromJSON(mash_details_path)
  mash_details_table <- data.frame(
    "Accession" = mash_details_data$reports$accession,
    "Identification" = mash_details_data$reports$organism$organism_name,
    "Tax ID" = mash_details_data$reports$organism$tax_id
  )

  # Top hits mash results
  mash_data <- read.csv(mash_path, sep="\t", header = F)
  
  pattern <- paste(mash_details_table$Accession, collapse = "|")
  mash_data_filtered <- mash_data %>%
    filter(str_detect(V1, pattern)) %>%
    mutate(V1 = str_extract(V1, pattern)) %>%
    rename("Accession" = V1,
           "Mash.distance" = V3,
           "Mash.p-value" = V4,
           "Mash.shared.hashes" = V5) %>%
    select(-V2)
  
  final_table <- mash_details_table %>%
    left_join(mash_data_filtered, by="Accession")
  
  final_table
}