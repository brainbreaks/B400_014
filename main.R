library(ggplot2)
library(cowplot)
library(ggseqlogo)
library(ShortRead)
library(Biostrings)


# pdf("htgts_logo.pdf", width=11.69, height=8.27, paper="a4r")
# sample_paths = list.files("/home/s215v/Workspace/B400_014_1_23447/preprocess/multx/", pattern="*.fq.gz", full.names=T)
# sample_paths = sample_paths[!grepl("unmatched", sample_paths)]
# plots = list()
# for(sample_path in sample_paths) {
#   sample_name = gsub("_.*", "", basename(sample_path))
#   pair_name = gsub(".*_(R[12]).*", "\\1", basename(sample_path))
#   r1 = ShortRead::readFastq(sample_path)
#   r1.reads = ShortRead::sread(r1)
#   r1.barcode = Biostrings::subseq(r1.reads, start=1, end=6)
#   r1.primer = Biostrings::subseq(r1.reads, start=7, end=30)
#   r1.adapter = Biostrings::subseq(r1.reads, start=31, end=49)
#   p1 = ggplot() +
#       ggseqlogo::geom_logo(as.character(r1.barcode)) +
#       labs(title="Barcode") +
#       ggseqlogo::theme_logo()
#   p2 = ggplot() +
#       ggseqlogo::geom_logo(as.character(r1.primer)) +
#       labs(title="Primer") +
#       ggseqlogo::theme_logo()
#   p3 = ggplot() +
#       ggseqlogo::geom_logo(as.character(r1.adapter)) +
#       labs(title="Adapter") +
#       ggseqlogo::theme_logo()
#   p = cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_label(paste0(sample_name, " (", pair_name, ")"), size=50), p1, p2, p3, ncol=1, title=sample_name)
#   plots[[sample_name]] = p
#   print(p)
# }
# dev.off()