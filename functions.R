
pacman::p_load('dplyr', 'tidyr', 'htmltools', 'bbplot', 'scales', 'data.table', 'ggmsa',
               'ggplot2', 'rdrop2', 'shiny', 'BiocManager', 'DECIPHER', 'ComplexUpset',
               'dendextend', 'data.table', 'Biostrings', 'alakazam', "unikn", 'ggupset',
               'plotly', "jcolors", 'ggdendro', "RColorBrewer","kmer","heatmaply", install = F)


load("2024-03-20_allele_igh_names.rda")
load("9_4_data_web.rda")

allele_appearance <- function(data_, imgt_genes, chain = "IGH") {
  
  data_ <- data_[grepl(imgt_genes, data_$gene),]
  
  alleles <- data_$allele
  data_$imgt_call <- data_$allele
  
  height = (length(unique(data_$imgt_call)))
  height = ifelse(length(height)<20, 25, height)
  height = height*30
  
  p <- ggplot(data_, aes(imgt_call)) + 
    geom_bar() +
    scale_x_discrete(drop=FALSE) +
    labs(x = "allele", y = "# Individuals", fill = "") + theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 10
      ),
      axis.text.y = element_text(
        size = 10
      ), axis.title = element_text(size = 14)
    )
  
  p1<-ggplotly(p, height = height, width = height)
  
  p <- ggplot(data_, aes(x = imgt_call, y = frac_allele,text = paste0(paste("</br>genotyped_allele : ",genotyped_allele,
                                                                            "</br>Subject: ", subject,
                                                                            "</br>lk: ", k_diff,
                                                                            "</br>in_genomic : ", in_genomic ,
                                                                            "</br>alleles: ",alleles,
                                                                            "</br>counts : ",counts ,
                                                                            "</br>genotyped_alleles : ",genotyped_alleles)))) +
    geom_boxplot() + # This adds the box plot # Flip coordinates if you still want horizontal boxes
    geom_point(aes(color = in_genomic), position = position_jitter(width = 0.2), alpha = 1, size = 1.5) +
    # Specify colors for 'yes' and 'no'
    scale_color_manual(values = c("yes" = alpha("darkblue", 1), "no" = alpha("#74A089", 1))) +
    labs(x = "Allele", y = "frec") + 
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12))
  
  p2<-ggplotly(p, height = height, width = height)
  
  sub_p1 <- subplot(p2, # boxplot and point - decrease size of points and remove legend
                    p1,
                    nrows = 2, margin = 0.007)
  
  return(sub_p1)
  
}




heatmap_alleles <-
  function(data_, g_group = "IGHVF5-G30", allele_db) {
    data_ <- data_[grepl(imgt_genes, data_$gene),]
    
    alleles <- data_$allele
    data_$imgt_call <- data_$allele
    data_upset <- data_ %>%
      select(subject, imgt_call, frac) %>%
      group_by(subject) %>%
      dplyr::mutate(frequency = sum(as.numeric(frac) / n(), na.rm = T),
                    val = 1, 
                    text = paste("</br>Subject: ", subject,
                                 "</br>Alleles: ", paste0(sort(unique(imgt_call)), collapse = ","),
                                 "</br>Frequency: ", unique(frequency)),
                    text2 = paste("</br>Alleles: ", paste0(sort(unique(imgt_call)), collapse = ","))) %>%
      select(-frac) %>%
      pivot_wider(names_from = imgt_call, 
                  values_from = val, 
                  values_fill = list(val = 0))
    
    
    
    t<-setNames(data_$imgt_call,data_$imgt_call)
    missing_alleles <- names(t)[!(t) %in% colnames(data_upset[,5:ncol(data_upset)])]
    
    data_upset[,missing_alleles] <- 0
    
    p_upset <- upset(
      data_upset,
      unique(t),
      annotations = list(
        'Frequency'=ggplot(mapping = aes(x=intersection, y = frequency)) +
          smplot2::sm_boxplot() + upset_themes$default
      ),
      width_ratio=0.1, encode_sets = F
    )
    
    p1 <- ggplotly(p_upset[[2]]+ aes(text = text), tooltip = "text")
    
    p2 <- ggplotly(p_upset[[4]] + aes(text = text2), tooltip = "text")
    p3 <- ggplotly(p_upset[[5]], tooltip = "text")
    p_interaction <- ggplotly(p_upset, tooltip = "x")
    
    sub_p1 <- subplot(hide_legend(p1), # boxplot and point - decrease size of points and remove legend
                      p2, 
                      hide_legend(p_interaction), 
                      nrows = 3, margin = 0.0007)
    sub_p2 <- subplot(plotly_empty(), 
                      plotly_empty(), 
                      p3, 
                      nrows = 3, margin = 0.0007)
    sub_p3 <- subplot(sub_p2, sub_p1, nrows = 1, widths= c(0.4,0.6) ,margin = 0.11)
    
    height = (length(unique(data_$imgt_call)))
    height = ifelse(length(height)<20, 25, height)
    height = height*100
    d<-ggplotly(sub_p3, height = height*1.5, width = height*1.2)
    return(d)
  }