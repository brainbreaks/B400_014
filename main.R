library(ggplot2)
library(cowplot)
library(ggseqlogo)
library(ShortRead)
library(Biostrings)
library(readr)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(GenomicRanges)
devtools::load_all('breaktools/')

effective_size = 1.87e9
sgRNA_length = 19
extsize = 1e5
maxgap = extsize*2
exttype = "symetrical"
threshold_qvalue = 1e-2
threshold_pileup = 1
slocal = 1e7
llocal = 1e7
bait_region=6e6
color_scheme = c("APH"="#1F78B4", "DMSO (APH)"="#A6CEE3")

samples_path = "~/Workspace/HTGTS_data/B400_014/tlx_samples.tsv"
samples_df = readr::read_tsv(samples_path) %>%
  dplyr::mutate(path=file.path(dirname(samples_path), path))

genes_ranges = rtracklayer::import("genomes/mm10/annotation/refGene.bed")
genes_ranges = genes_ranges[!grepl("_rev",genes_ranges$name)]
values(genes_ranges) = data.frame(gene_name=genes_ranges$name)
genes_reduced_ranges = GenomicRanges::reduce(genes_ranges)
genes_reduced_ranges$range_id = 1:length(genes_reduced_ranges)
genes_ranges = as.data.frame(IRanges::mergeByOverlaps(genes_ranges, genes_reduced_ranges)) %>%
  dplyr::arrange(dplyr::desc(genes_ranges.width)) %>%
  dplyr::distinct(range_id, .keep_all=T) %>%
  dplyr::mutate(seqnames=genes_reduced_ranges.seqnames, start=genes_reduced_ranges.start, end=genes_reduced_ranges.end) %>%
  dplyr::distinct(seqnames, start, end, gene_name) %>%
  GenomicRanges::makeGRangesFromDataFrame(ignore.strand=T, keep.extra.columns=T)


tlx_df = tlx_read_many(samples_df)
tlx_df = tlx_remove_rand_chromosomes(tlx_df)
tlx_df = tlx_mark_bait_chromosome(tlx_df)
tlx_df = tlx_mark_bait_junctions(tlx_df, bait_region)
baits_df = tlx_identify_baits(tlx_df, breaksite_size=sgRNA_length)
# tlx_df = tlx_mark_repeats(tlx_df, repeatmasker_df)
# tlx_df = tlx_df %>% dplyr::filter(!tlx_is_bait_junction)

sizes_df = readr::read_tsv("genomes/mm10/annotation/mm10.chrom.sizes", col_names=c("sizes_chrom", "sizes_length")) %>%
  dplyr::mutate(sizes_effective=sizes_length/sum(sizes_length)*effective_size)

libsizes_df = tlx_df %>%
  dplyr::group_by(tlx_group, tlx_sample, tlx_group_i, tlx_control) %>%
  dplyr::summarize(library_size=n()) %>%
  dplyr::mutate(Treatment=paste0(ifelse(tlx_control, "DMSO (", ""), tlx_group, ifelse(tlx_control, ")", ""))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(library_factor=max(library_size)/library_size)

ggplot(libsizes_df) +
  geom_bar(aes(x=tlx_group, y=library_size, fill=Treatment, group=paste0(tlx_control, tlx_group_i)), position="dodge", color="#EEEEEE", stat="identity") +
  geom_text(aes(x=tlx_group, y=library_size, group=paste0(tlx_control, tlx_group_i, Treatment), label=tlx_sample), vjust=-1, position=position_dodge(width=0.9), size=4) +
  scale_fill_manual(values=color_scheme) +
  scale_y_continuous(labels=scales::label_number(accuracy=1, scale=1e-3, suffix="K")) +
  labs(x="", y="Library size", title="Library size comparisons of Hydroxyurea/Aphidicolin treatment and their respective libraries") +
  theme_grey(base_size=14) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  theme(legend.position="bottom")


  #
  # 4. PCA samples
  #
  pca_binsize = 1e6
  tlx_hist_df = tlx_df %>%
    dplyr::filter(Rname=="chr6") %>%
    dplyr::inner_join(sizes_df, by=c("Rname"="sizes_chrom")) %>%
    dplyr::group_by(tlx_sample, Rname, sizes_length) %>%
    dplyr::do((function(d){
      dd<<-d
      h = hist(d$Junction, plot=F, breaks=c(seq(1, d$sizes_length[1], by=pca_binsize), d$sizes_length[1]))
      data.frame(tlx_sample=d$tlx_sample[1], hist_chrom=d$Rname[1], hist_start=h$breaks[-length(h$breaks)], hist_end=h$breaks[-1]-1, hist_count=h$count)
    })(.)) %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::arrange(dplyr::desc(hist_count)) %>%
    dplyr::mutate(is_top100=1:n()<=50) %>%
    dplyr::group_by(hist_chrom, hist_start, hist_end) %>%
    dplyr::filter(any(is_top100)) %>%
    dplyr::ungroup()


  tlx_hist_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_hist_df %>% dplyr::mutate(seqnames=Rname, start=hist_start, end=hist_end), ignore.strand=T, keep.extra.columns=T)
  tlx_named_hist_df = leftJoinByOverlaps(tlx_hist_ranges, genes_ranges) %>%
    dplyr::distinct(tlx_sample, hist_chrom, hist_start, hist_end, .keep_all=T) %>%
    dplyr::group_by(tlx_sample, hist_chrom, gene_name) %>%
    dplyr::summarize(hist_count=sum(hist_count))
  tlx_named_hist_wide_df = tlx_named_hist_df %>%
    dplyr::inner_join(libsizes_df %>% dplyr::select(tlx_sample, library_factor), by="tlx_sample") %>%
    dplyr::mutate(hist_count_norm=hist_count*library_factor) %>%
    reshape2::dcast(hist_chrom + gene_name ~ tlx_sample, value.var="hist_count_norm")
  tlx_mat = tlx_named_hist_wide_df %>%
    dplyr::mutate(rowname=paste(gene_name)) %>%
    tibble::column_to_rownames("rowname") %>%
    dplyr::select(-(hist_chrom:gene_name)) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("rowname") %>%
    dplyr::mutate(sample=rowname, rowname=sample) %>%
    dplyr::inner_join(samples_df, by="sample") %>%
    tibble::column_to_rownames("rowname") %>%
    dplyr::mutate(group_name=paste0(ifelse(control, "DMSO (", ""), group, ifelse(control, ")", "")))
  pca_res = prcomp(tlx_mat %>% dplyr::select(-(sample:group_name)), scale.=T)

  ggplot2::autoplot(pca_res, data=tlx_mat, colour="group_name", size=5, loadings=F, loadings.colour='blue', loadings.label=F, loadings.label.size=1) +
    ggrepel::geom_text_repel(aes(label=sample), size=6) +
    scale_color_manual(values=color_scheme) +
    labs(title=paste0("PCA of top 50 most abundant bins (", scales::number(pca_binsize, scale=1e-6,  suffix="Mbp"), ") from each sample"), color="Sample") +
    theme_grey(base_size=14)


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