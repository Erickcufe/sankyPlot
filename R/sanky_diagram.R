library(igraph)

snps <- vroom::vroom("data/pred_snps.csv")
snps <- snps[snps$value >= 2.33E-04,]

snps_info <- seekerBio::seeker_snp_context(snps$n2)
snps_arch <- seekerBio::seeker_snp_arq(snps$n2)

gd <- graph_from_data_frame(snps)

nodes <- get.vertex.attribute(gd)
nodes <- data.frame(nodes[[1]])
colnames(nodes) <- "name"
row.names(nodes)

nodes_select <- nodes$name
names(nodes_select) <- row.names(nodes)

links_1 <- data.frame(get.edgelist(gd))
colnames(links_1) <- c("source", "target")

for(i in 1:length(links_1$source)){
  b <- nodes_select[nodes_select == links_1$source[i]]
  links_1$source[i] <- names(b)
}

for(i in 1:length(links_1$target)){
  b <- nodes_select[nodes_select == links_1$target[i]]
  links_1$target[i] <- names(b)
}

links_1$source <- as.integer(links_1$source)
links_1$target <- as.integer(links_1$target)

links_2 <- cbind(links_1, value = edge.attributes(gd))
links_2$value <- -log2(links_2$value)
links_2 <- rbind(c(0, 1, 270), links_2)
links_2$source <- as.integer(links_2$source)
links_2$target <- as.integer(links_2$target)

gd_sanky <- list(nodes = nodes, links = links_2)


# Thus we can plot it
p <- sankeyNetwork(Links = gd_sanky$links, Nodes = gd_sanky$nodes, Source = "source")

p

# save the widget
library(htmlwidgets)
saveWidget(p, file=paste0( here::here(), "/sankeyClasification.html"))

# Load data

clasif <- vroom::vroom("data/clasificacion281masctl.csv")
risk <- vroom::vroom("data/risk_stimation_EUR.CSV")
bin1 <- vroom::vroom("data/Resultados_BIN1-2 (1).csv")


remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)


snps <- vroom::vroom("data/pred_snps.csv")
snps <- snps[snps$value >= 2.33E-04,]
snps$value <- NULL

context <- c("3_prime_UTR_variant", "regulatory_region_variant", "synonymous_variant")
genes <- c("CD36", "BIN1", "BIN1")
snps$context <- context
snps$genes <- genes
snps$BIN1 <- "BIN1"

readr::write_csv(snps, "sanky_data_snps.csv")

cell <- c("")

df <- readr::read_csv("sanky_data_snps.csv")
df <- df %>%
  make_long(n1, n2, context, genes, Cell, Expression)
df

library(ggplot2)
library(dplyr)

jpeg("sanky_snps.jpeg", units="in", width=10, height=10, res=300)
ggplot2::ggplot(df, ggplot2::aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = factor(node),
               label = node)) +
  geom_sankey(node.color = 1, flow.alpha = 0.7) +
  geom_sankey_label(color = 1, fill = "white") +
  theme_sankey(base_size = 16)+
  ggplot2::theme(legend.position = "none",
                 text = ggplot2::element_text(size = 25)) +
  ggplot2::scale_fill_viridis_d(option = "turbo")
dev.off()
