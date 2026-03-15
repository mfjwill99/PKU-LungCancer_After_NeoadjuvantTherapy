## R code for the FigureS1


################
##  load packages
################
library(GGally)
library(data.table)
library(ggsci)
library(gridExtra)
library(ggplot2)
library(ggsignif)
library(cowplot)

################
##  load data
################
##  FIGURE S1 A, B, C
T=read.csv("Data/T_Result_Index.csv")
##  FIGURE S1 D, E, F
LN=read.csv("Data/LN_Result_Index.csv")



################
##  load functions
################

#############################################
## function for distribution of LOH, TAI, LST
## Figure S1A, D
#############################################
my_custom_lower <- function(data, mapping, ...){
  ggplot(data = data, mapping = mapping) +
    geom_point(size = 1.5, alpha = 0.8) +#点的大小、颜色透明度
    #   geom_smooth(method = "lm", color = "red", linetype = "dashed", size=0.5) +#趋势线的颜色、类型、粗细
    theme_bw()
}

Distribution_plot<- function(df, 
                             outfile, 
                             plot_columns = c("LOH", "TAI", "LST"), 
                             group_col = "group",
                             colors = c("mLN+" = "#1F77B4", "mLN-" = "#FF7F0E", 
                                        "PT+" = "#FF7F0E", "PT-" = "#2CA02CFF"),
                             labels = c("mLN+" = "mLN+", "mLN-" = "mLN-", 
                                        "PT+" = "PT+", "PT-" = "PT-"),
                             lower_plot = NULL,
                             theme_settings = list(),
                             width = 6,
                             height = 5.5,
                             dpi = 300) {
  
  # Check required packages
  if (!requireNamespace("GGally", quietly = TRUE)) {
    stop("Package 'GGally' required. Install with: install.packages('GGally')")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required.")
  }
  
  # Validate inputs
  if (!all(plot_columns %in% colnames(df))) {
    stop("Some plot_columns not found in dataframe")
  }
  if (!group_col %in% colnames(df)) {
    stop(paste("Group column", group_col, "not found in dataframe"))
  }

# Default lower triangle plot function
if (is.null(lower_plot)) {
  my_custom_lower <- function(data, mapping, ...) {
    ggplot2::ggplot(data = data, mapping = mapping) + 
      ggplot2::geom_point(alpha = 0.7, size = 1) 
  }
  lower_plot <- my_custom_lower
}

# Merge user theme settings with defaults
default_theme <- list(
  strip.text.x = ggplot2::element_text(face = "bold", color = "black", size = 12, angle = 0),
  strip.text.y = ggplot2::element_text(face = "bold", color = "black", size = 12, angle = 0),
  axis.text.x = ggplot2::element_text(angle = 0),
  panel.grid = ggplot2::element_blank(),
  panel.background = ggplot2::element_rect(fill = "white"),
  plot.title = ggplot2::element_text(hjust = 0.5)
)
final_theme <- modifyList(default_theme, theme_settings)

# Create plot
p <- GGally::ggpairs(
  data = df,
  columns = plot_columns,
  mapping = ggplot2::aes_string(color = group_col),
  title = "",
  axisLabels = "show",
  lower = list(continuous = lower_plot),
  diag = list(continuous = "densityDiag", 
              mapping = ggplot2::aes_string(color = group_col)),
  upper = NULL
) + 
  ggplot2::scale_color_manual(values = colors, labels = labels) +
  ggplot2::scale_fill_manual(values = colors, labels = labels) +
  ggplot2::theme_bw() +
  do.call(ggplot2::theme, final_theme)

# Save plot
ggplot2::ggsave(
  filename = outfile,
  plot = p,
  width = width,
  height = height,
  dpi = dpi
)

message("Plot saved to: ", normalizePath(outfile))
invisible(p)
}


#############################################
## function for comparison
## Figure S1B,C,E,F
#############################################
compare_plot <- function(data, col = c("HRD", "LOH", "TAI", "LST"), 
                         colors = c("mLN+" = "#1F77B4", "mLN-" = "#FF7F0E", 
                                    "PT+" = "#FF7F0E", "PT-" = "#2CA02CFF"),
                         labels = c("mLN+" = "mLN+", "mLN-" = "mLN-", 
                                    "PT+" = "PT+", "PT-" = "PT-"),
                         outfile) {
  # 检查输入数据
  if (!"group" %in% colnames(data)) stop("Data must contain 'group' column")
  missing_cols <- setdiff(col, colnames(data))
  if (length(missing_cols) > 0) stop("Columns not found: ", paste(missing_cols, collapse = ", "))
  
  # 绘制单个变量的函数
  plot_single <- function(var) {
    df <- data.frame(group = data$group, value = data[[var]])
    df <- df[complete.cases(df), ]
    
    comparisons <- combn(unique(df$group), 2, simplify = FALSE)
    
    ggplot(df, aes(x = group, y = value, color = group)) +
      geom_boxplot(width = 0.5, outlier.shape = NA, size = 0.6) +
      geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
      stat_compare_means(
        comparisons = comparisons,
        method = "wilcox.test",
        label = "p.format",  # 改为显示具体p值
        tip.length = 0.01,
        size = 3,
        label.y = max(df$value, na.rm = TRUE) * 1.1  # 自动调整标注位置
      ) +
      scale_color_manual(values = colors) +  # 使用自定义颜色
      scale_fill_manual(values = colors) +   # 同步填充色（如有需要）
      labs(x = NULL, y = NULL, title = var) +
      theme_classic() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 9),
        axis.text.x = element_text(angle = 0, hjust = 1, size = 8),
        panel.border = element_rect(fill = NA, color = "black")
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))  # 为p值留出空间
  }
  
  # 生成所有子图
  plot_list <- lapply(col, plot_single)
  
  # 组合图形
  combined_plot <- cowplot::plot_grid(
    plotlist = plot_list,
    nrow = 1,
    align = "h",
    axis = "tb"
  )
  
  # 保存结果
  ggsave(
    outfile,
    combined_plot,
    width = length(col) * 2.5,
    height = 3.5,
    units = "in",
    device = "pdf"
  )
  
  message("Combined plot saved to: ", normalizePath(outfile))
  invisible(combined_plot)
}
################
##  Figure S1A
################
Distribution_plot(df = T, outfile = "FigureS1A_T.pdf",width = 6,height = 5.5)

################
##  Figure S1B
################
compare_plot(T, col = c("HRD", "LOH", "TAI", "LST"),  outfile = "FigureS1B_Compare.pdf")

################
##  Figure S1C
################
compare_plot(T, col = c("wgii.p.gain",	"wgii.t.gain", "wploidy"),  outfile = "FigureS1C_Compare.pdf")

################
##  Figure S1D
################
Distribution_plot(df = LN, outfile = "FigureS1D_LN.pdf",width = 6,height = 5.5)

################
##  Figure S1E
################
compare_plot(LN, col = c("HRD", "LOH", "TAI", "LST"),  outfile = "FigureS1E_Compare.pdf")

################
##  Figure S1F
################
compare_plot(LN, col = c("wgii.p.loss","wgii.t.loss"),  outfile = "FigureS1F_Compare.pdf")



