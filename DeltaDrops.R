# =====================================================================
# 愛媛県 言語変種データ解析スクリプト
# - δスコアに基づく系統的クラスタリングと地理空間可視化 -
# =====================================================================

# 1. パッケージの読み込み ----------------------------------------------------
library(tidyverse)
library(sf)
library(deldir)
library(igraph)
library(geosphere)
library(ggrepel)
library(mapdata)

# 2. 設定 -------------------------------------------------------------------
# 作業ディレクトリとファイル名を指定してください。
target_dir <- "./data/" # 入力データの場所
target_dir <- "./output/" # 出力ファイルの保存先
input_file <- "" # 入力データファイル名を書いてください．
output_prefix <- ""　# 出力ファイルのプレフィックスを書いてください．

# パラメータ設定
dist_buffer_km <- 30 # Delaunay 網のエッジを切断する距離（km）
delta_threshold <- 0.486 # δスコアの閾値
voronoi_buffer_km <- 20 # Voronoi 分割の境界クリップ用バッファ（km）

# 3. データの読み込みと前処理 -------------------------------------------------
file_path <- file.path(target_dir, input_file)

# 読み込み（ファイルが存在しない場合はエラーを出して停止します）
if(!file.exists(file_path)) {
  stop(paste("エラー: データファイルが見つかりません ->", file_path))
}

data <- read_tsv(file_path, show_col_types = FALSE)

# Taxon, Lat, Longの分離と形質行列の作成
taxa_names <- data$Taxon
coords <- data %>% select(Long, Lat)
trait_matrix <- data %>% select(-Taxon, -Lat, -Long) %>% as.matrix()
rownames(trait_matrix) <- taxa_names

# 4. 形質距離行列（0,1なのでマンハッタン距離でやっています）の作成 -----------------------------------
dist_trait <- dist(trait_matrix, method = "manhattan")
dist_trait_mat <- as.matrix(dist_trait)

# 5. 地理空間ネットワーク（Delaunay網 + 30kmバッファ）の作成 ------------------
# deldirでDelaunay網を作成
del <- deldir(coords$Long, coords$Lat)
edges <- del$delsgs %>% select(ind1, ind2)

# 各エッジの地理的距離を計算
edges$dist_m <- apply(edges, 1, function(row) {
  p1 <- c(coords$Long[row[1]], coords$Lat[row[1]])
  p2 <- c(coords$Long[row[2]], coords$Lat[row[2]])
  distGeo(p1, p2)
})

# 30km (30,000m) 以内のエッジのみを抽出し、隣接グラフを作成
edges_filtered <- edges %>% filter(dist_m <= dist_buffer_km * 1000)
adj_graph <- graph_from_data_frame(
  d = data.frame(
    from = taxa_names[edges_filtered$ind1], 
    to = taxa_names[edges_filtered$ind2]
  ), 
  directed = FALSE, 
  vertices = taxa_names
)

# 6. δスコア計算関数の定義（スクラッチ実装） ----------------------------------
# 4点からδスコアを計算する関数
calc_quartet_delta <- function(d_mat, t1, t2, t3, t4) {
  s1 <- d_mat[t1, t2] + d_mat[t3, t4]
  s2 <- d_mat[t1, t3] + d_mat[t2, t4]
  s3 <- d_mat[t1, t4] + d_mat[t2, t3]
  
  sums <- sort(c(s1, s2, s3), decreasing = TRUE)
  L1 <- sums[1]
  L2 <- sums[2]
  L3 <- sums[3]
  
  if ((L1 - L3) == 0) {
    return(0)
  } else {
    return((L1 - L2) / (L1 - L3))
  }
}

# 任意の数の変種集合からδスコア（全4点の平均）を計算する関数
calc_group_delta <- function(d_mat, taxa_set) {
  if (length(taxa_set) < 4) return(NA)
  
  quartets <- combn(taxa_set, 4)
  deltas <- apply(quartets, 2, function(q) {
    calc_quartet_delta(d_mat, q[1], q[2], q[3], q[4])
  })
  
  return(mean(deltas))
}

# 7. 初期コアの選定 ----------------------------------------------------------
# Delaunay網上で連結しているすべての4点部分グラフを列挙
all_quartets <- combn(taxa_names, 4)
connected_quartets_idx <- apply(all_quartets, 2, function(q) {
  sub_g <- induced_subgraph(adj_graph, q)
  is_connected(sub_g)
})

valid_quartets <- all_quartets[, connected_quartets_idx]

# 各4点のδスコアを計算し、評価
core_eval <- data.frame(
  id = 1:ncol(valid_quartets),
  taxa = apply(valid_quartets, 2, paste, collapse = ","),
  delta = apply(valid_quartets, 2, function(q) calc_group_delta(dist_trait_mat, q))
)

# 閾値未満で、スコアが低い順にソート
core_candidates <- core_eval %>%
  filter(delta < delta_threshold) %>%
  arrange(delta) %>%
  mutate(taxa_list = strsplit(as.character(taxa), ","))

# 8. グローバル貪欲法による拡張アルゴリズム ----------------------------------
expansion_history <- data.frame(
  Group_ID = integer(),
  Size = integer(),
  Delta_Score = numeric(),
  Added_Taxon = character(),
  stringsAsFactors = FALSE
)

final_groups <- list()
used_cores <- c()
group_counter <- 1

# 拡張処理のメインループ
for (i in 1:nrow(core_candidates)) {
  current_core <- core_candidates$taxa_list[[i]]
  
  # すでに確立した別の方言群と全く重複していないかチェック（完全独立のコアのみ採用）
  is_overlap <- any(sapply(final_groups, function(g) length(intersect(g, current_core)) > 0))
  if (is_overlap) next
  
  # 新しい方言群として拡張開始
  current_group <- current_core
  
  # 初期コアの履歴を記録（4つまとめて最初のステップとする）
  expansion_history <- rbind(expansion_history, data.frame(
    Group_ID = group_counter,
    Size = 4,
    Delta_Score = core_candidates$delta[i],
    Added_Taxon = paste(current_group, collapse = ",")
  ))
  
  # 拡張ループ
  repeat {
    # 現在の群に隣接し、かつ現在の群にまだ含まれていない変種のリストを取得
    # ※他の群に属しているか（オーバーラップ）はここでは問わず吸収テストを行う
    neighbors <- unlist(lapply(current_group, function(v) neighbors(adj_graph, v)$name))
    candidate_taxa <- unique(neighbors[!neighbors %in% current_group])
    
    if (length(candidate_taxa) == 0) break # 拡張可能な隣接変種がない
    
    best_candidate <- NULL
    min_delta <- Inf
    
    # 各候補を仮に追加してδスコアを計算
    for (cand in candidate_taxa) {
      test_group <- c(current_group, cand)
      test_delta <- calc_group_delta(dist_trait_mat, test_group)
      
      if (test_delta < delta_threshold && test_delta < min_delta) {
        min_delta <- test_delta
        best_candidate <- cand
      }
    }
    
    # 閾値未満で追加できる候補がいれば吸収、いなければ終了
    if (!is.null(best_candidate)) {
      current_group <- c(current_group, best_candidate)
      expansion_history <- rbind(expansion_history, data.frame(
        Group_ID = group_counter,
        Size = length(current_group),
        Delta_Score = min_delta,
        Added_Taxon = best_candidate
      ))
    } else {
      break
    }
  }
  
  final_groups[[group_counter]] <- current_group
  group_counter <- group_counter + 1
}

# 9. 変種の所属カテゴリ分類 --------------------------------------------------
taxon_membership <- data.frame(Taxon = taxa_names, stringsAsFactors = FALSE)
taxon_membership$groups <- lapply(taxa_names, function(t) {
  which(sapply(final_groups, function(g) t %in% g))
})

taxon_membership <- taxon_membership %>%
  mutate(
    group_count = sapply(groups, length),
    category = case_when(
      group_count == 1 ~ "Single Group",
      group_count > 1 ~ "Boundary Zone (Overlap)",
      group_count == 0 ~ "Isolated"
    ),
    primary_group = ifelse(group_count > 0, sapply(groups, function(g) g[1]), NA)
  )

# 地図用のデータとマージ
map_data_df <- left_join(data, taxon_membership, by = "Taxon")

# 10. ボロノイ地図の作成 -----------------------------------------------------
# sfオブジェクトに変換（日本の平面直角座標系 EPSG:6672 等を用いると正確ですが、
# 汎用性を持たせるため地理座標系 EPSG:4326 から UTM (EPSG:32653) を経由します）
sf_points <- st_as_sf(map_data_df, coords = c("Long", "Lat"), crs = 4326) %>%
  st_transform(crs = 32653)

# バッファで囲む領域を作成
boundary_buffer <- st_buffer(st_union(sf_points), dist = voronoi_buffer_km * 1000)

# ボロノイ分割を作成し、バッファでクリップ
voronoi_polys <- st_voronoi(st_union(sf_points)) %>%
  st_collection_extract("POLYGON") %>%
  st_sf() %>%
  st_join(sf_points) %>%
  st_intersection(boundary_buffer) %>%
  st_transform(crs = 4326)

sf_points <- st_transform(sf_points, crs = 4326)

# 色の設定（方言群には青系、境界には赤系、孤立にはグレー）
num_groups <- max(0, length(final_groups))
blue_palette <- colorRampPalette(c("#E0F3F8", "#4575B4", "#313695"))(num_groups)

# マッピング用の色付け列を作成
voronoi_polys <- voronoi_polys %>%
  mutate(fill_color = case_when(
    category == "Isolated" ~ "#CCCCCC",
    category == "Boundary Zone (Overlap)" ~ "#D73027",
    category == "Single Group" ~ blue_palette[primary_group]
  ))

# 海岸線データ
japan_coast <- map_data("japan")

p_map <- ggplot() +
  geom_polygon(data = japan_coast, aes(x = long, y = lat, group = group), 
               fill = "#F0F0F0", color = "gray50", size = 0.2) +
  geom_sf(data = voronoi_polys, aes(fill = fill_color), alpha = 0.7, color = "white") +
  geom_sf(data = sf_points, size = 1, color = "black") +
  geom_text_repel(data = map_data_df, aes(x = Long, y = Lat, label = Taxon),
                  size = 3, bg.color = "white", bg.r = 0.15) +
  scale_fill_identity() +
  coord_sf(xlim = range(map_data_df$Long) + c(-0.2, 0.2), 
           ylim = range(map_data_df$Lat) + c(-0.2, 0.2), expand = FALSE) +
  theme_minimal() +
  labs(title = "Language Varieties: Voronoi Map of Dialect Groups",
       x = "Longitude", y = "Latitude")

# 地図の保存
ggsave(file.path(target_dir, paste0(output_prefix, "_Voronoi_map.png")), 
       plot = p_map, width = 10, height = 8, dpi = 300)

# 11. δスコア推移の折れ線グラフ作成 ------------------------------------------
# 折れ線グラフ用にデータを加工
plot_data <- expansion_history %>%
  mutate(Group_Factor = as.factor(Group_ID))

p_traj <- ggplot(plot_data, aes(x = Size, y = Delta_Score, color = Group_Factor, group = Group_Factor)) +
  geom_hline(yintercept = delta_threshold, linetype = "dashed", color = "red", size = 0.8) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = Added_Taxon), size = 3, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = blue_palette) +
  theme_minimal() +
  labs(title = "Delta Score Trajectory During Group Expansion",
       x = "Number of Varieties in Group",
       y = expression(delta ~ "Score"),
       color = "Group ID")

# グラフの保存
ggsave(file.path(target_dir, paste0(output_prefix, "_trajectory_plot.png")), 
       plot = p_traj, width = 10, height = 6, dpi = 300)

# 12. 拡張履歴のCSV出力 ------------------------------------------------------
write_csv(expansion_history, file.path(target_dir, paste0(output_prefix, "_history.csv")))

# 終了メッセージ
message("解析が完了しました。以下のファイルが出力されました：")
message("1. ", file.path(target_dir, paste0(output_prefix, "_Voronoi_map.png")))
message("2. ", file.path(target_dir, paste0(output_prefix, "_trajectory_plot.png")))
message("3. ", file.path(target_dir, paste0(output_prefix, "_history.csv")))