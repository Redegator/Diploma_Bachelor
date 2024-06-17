setwd("D:/Diploma")  # Встановлення робочої директорії
getwd()  # Виведення поточної робочої директорії

# Завантаження необхідних бібліотек
library(readxl)  
library(openxlsx)  
library(dplyr)  
library(AnnotationDbi)  
library(clusterProfiler)
library(enrichplot)  
library(org.Hs.eg.db)
library(ggplot2)

# Завантаження файлу з матрицею експресії
count_matrix <- read.table("df_merged.txt", header = TRUE, row.names = 1)

# Видалення останніх 5 рядків з count_matrix
count_matrix <- head(count_matrix, -5)

# Вказання шляху до файлу Excel
file_path <- "suzuki_data.xlsx"

# Зчитування файлу Excel у датафрейм
results <- read_excel(file_path)

colnames(results)  # Виведення назв колонок

# Зміна знаку значень в колонці з logFC
results$M <- results$M * -1
results[c('Official gene symbol','M')]

# Вибір генів з позитивною регуляцією
up_reg <- results[which(results$M > 1),]
# Вибір генів з негативною регуляцією
down_reg <- results[which(results$M < -1),]

# Сортування генів з позитивною регуляцією за спаданням logFC
up_reg <- up_reg[order(up_reg$M, decreasing = TRUE), ]
# Сортування генів з негативною регуляцією за спаданням logFC
down_reg <- down_reg[order(down_reg$M, decreasing = TRUE),]

# Створення нового файлу Excel
wb <- createWorkbook()

# Додавання нового листа у файл
addWorksheet(wb, "Data")
# Запис даних у лист
writeData(wb, "Data", up_reg[c("Ensemble IDs", 'Official gene symbol','M')])

# Збереження файлу Excel
saveWorkbook(wb, "Suzuki_DEGs.xlsx", overwrite = TRUE)

# Вибір Ensemble ID генів
DEGs <- results$"Ensemble IDs"
# Вибір всіх імен генів з матриці експресії
all_genes <- rownames(count_matrix)
length(all_genes)  # Виведення кількості всіх генів

# Створення логічного вектора, який вказує, чи всі значення в рядку не дорівнюють нулю
non_zero_genes <- apply(count_matrix, 1, function(x) all(x != 0))

# Сортування результатів за значеннями M
sorted_results <- arrange(results, M)
sorted_results <- arrange(results, desc(M))

# Вибір генів з позитивною регуляцією (logFC > 2)
up_reg <- results[which(results$M  > 2),]
# Вибір генів з негативною регуляцією (logFC < -2)
down_reg <- results[which(results$M < -2),]

organism <- "org.Hs.eg.db"  # Визначення організму для бази даних анотацій

# Отримання назв генів і типів на основі ENSG ID
gene_info <- AnnotationDbi::select(org.Hs.eg.db, keys = DEGs, columns = c("SYMBOL", "GENETYPE"), keytype = "ENSEMBL")
gene_info$logFC <- results$M  # Додавання стовпця logFC до gene_info

# Об'єднання двох датафреймів
combined_data <- rbind(up_reg, down_reg)

# Аналіз збагачення GO
go_enrich <- enrichGO(gene = DEGs,
                      universe = all_genes,
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)



# Ініціалізація змінної для підрахунку ДЕГів, які не повертають результатів
c = 0
# Ініціалізація вектора для збереження індексів цих ДЕГів
error_indices <- c()
# Перебір кожного гена з DEGs
for (DEG in DEGs) {
  c = c + 1
  print(c)
  
  # Аналіз збагачення GO для кожного гена
  go_e <- enrichGO(gene = DEG,
                   universe = all_genes,
                   OrgDb = organism, 
                   keyType = 'ENSEMBL',
                   readable = T,
                   ont = "BP",
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.10)
  
  # Якщо результату немає, індексу ДЕГа зберігається
  if (is.null(go_e)) {
    error_indices <- append(error_indices, c)
  }
}

# Збереження індексів ДЕГів, які не повертають результатів
error_indices <- c(6, 20, 21, 23, 24, 27, 28, 33, 34, 35, 36, 37, 38, 39, 40, 41)
DEGs[error_indices]  # Виведення генів з помилками

# Виведення інформації про ДЕГи, які не повертають результатів
gene_info[error_indices, ]

# Всі індекси датафрейму
all_indices <- 1:nrow(gene_info)

# Індекси, які не містяться в error_indices
valid_indices <- setdiff(all_indices, error_indices)

# Виведення рядків з позицій, яких нема у векторі error_indices
gene_info[valid_indices, ]

length(error_indices)  # Виведення кількості ДЕГів, які не повертають результатів
length(valid_indices)  # Виведення кількості ДЕГів, які не повертають результатів

# Додавання специфічних індексів до валідних індексів
valid_indices <- c(valid_indices, 27, 28)

# Сортування інформації про гени за значенням logFC
sorted_gene_info <- gene_info[valid_indices, ][order(gene_info[valid_indices, ]$logFC, 
                                                     decreasing = TRUE),]

length(rownames(sorted_gene_info))  # Виведення кількості рядків у sorted_gene_info

# Встановлення кольорів для регуляції
colors <- ifelse(sorted_gene_info$logFC > 0, 'green', 'red')

par(mar = c(6, 5, 2, 2) + 0.1)  # Налаштування полів для графіку

# Побудова стовпчикової діаграми ДЕГів і їх logFC
barplot(height = sorted_gene_info$logFC, 
        names.arg = sorted_gene_info$SYMBO, 
        col = colors, 
        ylab = "log2FC",
        las = 2)

# Додавання легенди до графіку
legend("topright", legend = c("Позитивно регульовані", "Негативно регульовані"), fill = c("green", "red"))

# Вибір тестових ДЕГів
DEGs_test <- DEGs[c(valid_indices, 27, 28)]

# Перевірковий аналіз збагачення GO для тестових ДЕГів
go_enrich_test_half <- enrichGO(gene = DEGs_test,
                                universe = all_genes,
                                OrgDb = organism, 
                                keyType = 'ENSEMBL',
                                readable = T,
                                ont = "BP",
                                pvalueCutoff = 0.05, 
                                qvalueCutoff = 0.10)

# Виведення інформації про гени з валідних індексів та специфічних індексів
gene_info[c(valid_indices, 27, 28), ]

# Запис даних у лист
writeData(wb, "Data", cluster_summary[c("Description", 'Count','p.adjust', 'geneID')])

# Збереження файлу Excel
saveWorkbook(wb, "Suzuki_pathways.xlsx", overwrite = TRUE)

# Побудова графіку dot plot для GO збагачення
dotplot(go_enrich, showCategory = 19)

# Обчислення матриці термінової подібності
pairwise <- pairwise_termsim(go_enrich)
# Побудова графіку emap plot
emapplot(pairwise)

# Встановлення значень logFC для кожного гена
foldChange <- setNames(results$M, results$`Ensemble IDs`)

# Побудова графіку cnet plot
cnetplot(go_enrich, showCategory = 8,
         color.params = list(foldChange = foldChange))

# Побудова графіку з кольорами на основі foldChange
p <- cnetplot(go_enrich, showCategory = 8, 
              color.params = list(foldChange = foldChange),
              cex_label_category = 0.8)

# Отримання даних з графіку
p_data <- p$data

# Ідентифікація вузлів, що є шляхами (категоріями)
pathway_nodes <- p_data[p_data$name %in% go_enrich$Description, ]

# Оновлення кольору вузлів шляхів
p <- p + 
  geom_point(data = pathway_nodes, aes(x = x, y = y, size = size), color = "#e57373")

# Додавання підписів до графіку
p <- p + 
  labs(size = "Кількість генів у процесі", color = "log2FC")

# Додавання легенди до графіку
p <- p +
  guides(color = guide_colorbar(title = "log2FC"))

# Відображення оновленого графіку
print(p)
