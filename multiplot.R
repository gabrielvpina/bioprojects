# Carregar pacotes
library(ggbio)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(cowplot)  # Para arranjar os múltiplos gráficos

# Carregar arquivos GFF3
gff_files <- list.files(path = "~/laiana/annot_normalized", pattern = "*.gff3", full.names = TRUE)

# Função para carregar e processar os arquivos GFF3
gff_data <- lapply(gff_files, function(file) {
  gff <- import(file)  # Importa o arquivo GFF3
  return(gff)
})

# Filtrar as regiões CDS e genes de cada arquivo
cds_genes <- lapply(gff_data, function(gff) {
  cds <- gff[which(gff$type %in% c("gene", "CDS")), ]
  return(cds)
})

# Função para criar o gráfico de genes/CDS para um GFF3
create_plot <- function(cds) {
  # Criar o gráfico utilizando ggplot
  p <- ggplot() +
    geom_rect(data = as.data.frame(cds), 
              aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = type)) +
    scale_fill_manual(values = c("gene" = "blue", "CDS" = "red")) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x = "Posição", y = "Genoma") +
    scale_x_continuous(expand = c(0, 0)) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  return(p)
}

# Criar os gráficos para todos os arquivos GFF3
plots <- lapply(cds_genes, create_plot)

# Combinar todos os gráficos em um único painel
combined_plot <- plot_grid(plotlist = plots, ncol = 1)

# Exibir o gráfico combinado
print(combined_plot)

# Salvar o gráfico combinado em PDF
ggsave("genes_cds_combined_plot.pdf", plot = combined_plot, width = 10, height = length(gff_files) * 3)


