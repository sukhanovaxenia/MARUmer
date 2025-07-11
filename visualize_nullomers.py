#!/usr/bin/env python3
"""
visualize_nullomers.py - Визуализация и углубленный анализ нулломеров
"""

import argparse
import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import networkx as nx
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Настройка стиля графиков
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

def load_debruijn_analysis(filepath):
    """Загрузка результатов анализа графа де Брюйна"""
    nullomers = []
    
    with open(filepath, 'r') as f:
        in_details = False
        for line in f:
            if line.startswith('nullomer\tclassification'):
                in_details = True
                continue
            
            if in_details and line.strip() and not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) >= 7:
                    nullomer_data = {
                        'nullomer': parts[0],
                        'classification': parts[1],
                        'gc_content': float(parts[2]),
                        'pred_count': int(parts[3]),
                        'succ_count': int(parts[4]),
                        'predecessors': parts[5].split(',') if parts[5] != 'none' else [],
                        'successors': parts[6].split(',') if parts[6] != 'none' else []
                    }
                    nullomers.append(nullomer_data)
    
    return pd.DataFrame(nullomers)

def analyze_motifs(df, k=13):
    """Анализ мотивов в нулломерах"""
    print("Анализ мотивов в нулломерах...")
    
    # Поиск перепредставленных подпоследовательностей
    motif_lengths = [3, 4, 5]
    motif_stats = {}
    
    for length in motif_lengths:
        motif_counts = defaultdict(int)
        
        for nullomer in df['nullomer']:
            for i in range(len(nullomer) - length + 1):
                motif = nullomer[i:i+length]
                motif_counts[motif] += 1
        
        # Находим топ-10 самых частых мотивов
        top_motifs = sorted(motif_counts.items(), key=lambda x: x[1], reverse=True)[:10]
        motif_stats[length] = top_motifs
    
    return motif_stats

def create_comprehensive_visualizations(df, output_dir, k=13):
    """Создание комплексного набора визуализаций"""
    
    # 1. Расширенный анализ классификаций
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Круговая диаграмма классификаций
    class_counts = df['classification'].value_counts()
    axes[0, 0].pie(class_counts.values, labels=class_counts.index, autopct='%1.1f%%')
    axes[0, 0].set_title('Распределение нулломеров по классификациям')
    
    # GC-содержание по классификациям
    df.boxplot(column='gc_content', by='classification', ax=axes[0, 1])
    axes[0, 1].set_title('GC-содержание по типам нулломеров')
    axes[0, 1].set_xlabel('Классификация')
    axes[0, 1].set_ylabel('GC-содержание')
    
    # Связность нулломеров
    df['total_connections'] = df['pred_count'] + df['succ_count']
    axes[1, 0].hist(df['total_connections'], bins=20, edgecolor='black')
    axes[1, 0].set_xlabel('Общее количество связей')
    axes[1, 0].set_ylabel('Количество нулломеров')
    axes[1, 0].set_title('Распределение связности нулломеров')
    
    # Корреляция между предшественниками и последователями
    axes[1, 1].scatter(df['pred_count'], df['succ_count'], alpha=0.5)
    axes[1, 1].set_xlabel('Количество предшественников')
    axes[1, 1].set_ylabel('Количество последователей')
    axes[1, 1].set_title('Корреляция между входящими и исходящими связями')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/comprehensive_analysis.png', dpi=300)
    plt.close()
    
    # 2. Анализ длин повторов
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Подсчет гомополимеров
    homopolymer_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    max_repeat_lengths = []
    
    for nullomer in df['nullomer']:
        max_repeat = 1
        current_repeat = 1
        
        for i in range(1, len(nullomer)):
            if nullomer[i] == nullomer[i-1]:
                current_repeat += 1
                max_repeat = max(max_repeat, current_repeat)
            else:
                current_repeat = 1
        
        max_repeat_lengths.append(max_repeat)
        
        # Подсчет гомополимеров длины 3+
        for base in 'ACGT':
            if base * 3 in nullomer:
                homopolymer_counts[base] += 1
    
    # График максимальных длин повторов
    repeat_hist = np.bincount(max_repeat_lengths)
    ax.bar(range(len(repeat_hist)), repeat_hist)
    ax.set_xlabel('Максимальная длина повтора')
    ax.set_ylabel('Количество нулломеров')
    ax.set_title('Распределение максимальных длин повторов в нулломерах')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/repeat_analysis.png', dpi=300)
    plt.close()
    
    # 3. Тепловая карта динуклеотидных частот
    dinuc_matrix = np.zeros((4, 4))
    bases = 'ACGT'
    base_to_idx = {b: i for i, b in enumerate(bases)}
    
    for nullomer in df['nullomer']:
        for i in range(len(nullomer) - 1):
            if nullomer[i] in bases and nullomer[i+1] in bases:
                dinuc_matrix[base_to_idx[nullomer[i]], base_to_idx[nullomer[i+1]]] += 1
    
    # Нормализация
    dinuc_matrix = dinuc_matrix / dinuc_matrix.sum()
    
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(dinuc_matrix, annot=True, fmt='.3f', cmap='YlOrRd',
                xticklabels=list(bases), yticklabels=list(bases))
    ax.set_xlabel('Вторая база')
    ax.set_ylabel('Первая база')
    ax.set_title('Частоты динуклеотидов в нулломерах')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/dinucleotide_frequencies.png', dpi=300)
    plt.close()
    
    # 4. PCA анализ k-мерного состава
    print("Выполнение PCA анализа...")
    
    # Создание матрицы признаков
    feature_matrix = []
    for nullomer in df['nullomer']:
        features = []
        # Состав оснований
        for base in 'ACGT':
            features.append(nullomer.count(base) / k)
        # Динуклеотидный состав
        for i in range(4):
            for j in range(4):
                dinuc = bases[i] + bases[j]
                features.append(nullomer.count(dinuc) / (k-1))
        feature_matrix.append(features)
    
    X = np.array(feature_matrix)
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X)
    
    # Визуализация PCA с раскраской по классификации
    fig, ax = plt.subplots(figsize=(10, 8))
    
    for classification in df['classification'].unique():
        mask = df['classification'] == classification
        ax.scatter(X_pca[mask, 0], X_pca[mask, 1], label=classification, alpha=0.6)
    
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    ax.set_title('PCA проекция нулломеров')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/pca_analysis.png', dpi=300)
    plt.close()
    
    # 5. Анализ мотивов
    motif_stats = analyze_motifs(df, k)
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    for idx, (length, motifs) in enumerate(motif_stats.items()):
        motif_names = [m[0] for m in motifs[:10]]
        motif_counts = [m[1] for m in motifs[:10]]
        
        axes[idx].barh(motif_names, motif_counts)
        axes[idx].set_xlabel('Частота')
        axes[idx].set_title(f'Топ-10 {length}-меров в нулломерах')
        axes[idx].invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/motif_analysis.png', dpi=300)
    plt.close()

def create_network_visualization(df, output_dir, sample_size=100):
    """Создание визуализации локальной структуры графа вокруг нулломеров"""
    print("Создание сетевой визуализации...")
    
    # Выбираем случайную выборку нулломеров для визуализации
    sample_df = df.sample(min(sample_size, len(df)))
    
    G = nx.DiGraph()
    
    # Добавляем узлы и ребра
    for _, row in sample_df.iterrows():
        nullomer = row['nullomer']
        
        # Добавляем нулломер как отсутствующий узел
        G.add_node(nullomer, node_type='nullomer', exists=False)
        
        # Добавляем предшественников
        for pred in row['predecessors']:
            G.add_node(pred, node_type='kmer', exists=True)
            G.add_edge(pred, nullomer, edge_type='to_nullomer')
        
        # Добавляем последователей
        for succ in row['successors']:
            G.add_node(succ, node_type='kmer', exists=True)
            G.add_edge(nullomer, succ, edge_type='from_nullomer')
    
    # Визуализация
    plt.figure(figsize=(20, 20))
    
    pos = nx.spring_layout(G, k=2, iterations=50)
    
    # Раскраска узлов
    node_colors = []
    for node in G.nodes():
        if G.nodes[node]['node_type'] == 'nullomer':
            node_colors.append('red')
        else:
            node_colors.append('lightblue')
    
    # Рисуем граф
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=300, alpha=0.8)
    nx.draw_networkx_edges(G, pos, edge_color='gray', arrows=True, alpha=0.5)
    
    plt.title(f'Локальная структура графа де Брюйна вокруг {sample_size} нулломеров')
    plt.axis('off')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/network_structure.png', dpi=300, bbox_inches='tight')
    plt.close()

def generate_summary_report(df, output_dir, k=13):
    """Генерация сводного отчета"""
    report_path = f'{output_dir}/summary_report.txt'
    
    with open(report_path, 'w') as f:
        f.write("=== СВОДНЫЙ ОТЧЕТ ПО АНАЛИЗУ НУЛЛОМЕРОВ ===\n\n")
        
        # Основная статистика
        f.write(f"Размер k-мера: {k}\n")
        f.write(f"Всего нулломеров проанализировано: {len(df)}\n")
        f.write(f"Всего возможных {k}-меров: {4**k:,}\n")
        f.write(f"Процент нулломеров: {100 * len(df) / (4**k):.4f}%\n\n")
        
        # Статистика по классификациям
        f.write("КЛАССИФИКАЦИЯ НУЛЛОМЕРОВ:\n")
        class_stats = df['classification'].value_counts()
        for classification, count in class_stats.items():
            f.write(f"  {classification}: {count} ({100 * count / len(df):.1f}%)\n")
        f.write("\n")
        
        # GC-содержание
        f.write("GC-СОДЕРЖАНИЕ:\n")
        f.write(f"  Среднее: {df['gc_content'].mean():.3f}\n")
        f.write(f"  Медиана: {df['gc_content'].median():.3f}\n")
        f.write(f"  Станд. откл.: {df['gc_content'].std():.3f}\n")
        f.write(f"  Минимум: {df['gc_content'].min():.3f}\n")
        f.write(f"  Максимум: {df['gc_content'].max():.3f}\n\n")
        
        # Связность
        f.write("СВЯЗНОСТЬ В ГРАФЕ:\n")
        f.write(f"  Среднее кол-во предшественников: {df['pred_count'].mean():.2f}\n")
        f.write(f"  Среднее кол-во последователей: {df['succ_count'].mean():.2f}\n")
        f.write(f"  Изолированных нулломеров: {len(df[df['classification'] == 'isolated'])}\n")
        f.write(f"  Нулломеров без предшественников: {len(df[df['pred_count'] == 0])}\n")
        f.write(f"  Нулломеров без последователей: {len(df[df['succ_count'] == 0])}\n\n")
        
        # Экстремальные случаи
        f.write("ЭКСТРЕМАЛЬНЫЕ СЛУЧАИ:\n")
        
        # Самые связанные нулломеры
        df['total_connections'] = df['pred_count'] + df['succ_count']
        most_connected = df.nlargest(5, 'total_connections')
        f.write("  Наиболее связанные нулломеры:\n")
        for _, row in most_connected.iterrows():
            f.write(f"    {row['nullomer']}: {row['total_connections']} связей\n")
        
        # Экстремальные GC
        f.write("\n  Экстремальное GC-содержание:\n")
        lowest_gc = df.nsmallest(3, 'gc_content')
        highest_gc = df.nlargest(3, 'gc_content')
        
        f.write("    Минимальное GC:\n")
        for _, row in lowest_gc.iterrows():
            f.write(f"      {row['nullomer']}: {row['gc_content']:.3f}\n")
        
        f.write("    Максимальное GC:\n")
        for _, row in highest_gc.iterrows():
            f.write(f"      {row['nullomer']}: {row['gc_content']:.3f}\n")
    
    print(f"Сводный отчет сохранен в {report_path}")

def main():
    parser = argparse.ArgumentParser(description='Визуализация и анализ результатов поиска нулломеров')
    parser.add_argument('results_dir', help='Директория с результатами анализа')
    parser.add_argument('--k', type=int, default=13, help='Размер k-мера')
    parser.add_argument('--format', choices=['debruijn', 'jellyfish'], default='debruijn',
                       help='Формат входных данных')
    
    args = parser.parse_args()
    
    # Поиск файла с результатами
    if args.format == 'debruijn':
        result_file = os.path.join(args.results_dir, 'debruijn_analysis.txt')
    else:
        result_file = os.path.join(args.results_dir, 'nullomer_analysis.txt')
    
    if not os.path.exists(result_file):
        print(f"Файл {result_file} не найден!")
        return
    
    # Загрузка данных
    print(f"Загрузка данных из {result_file}...")
    df = load_debruijn_analysis(result_file)
    
    if df.empty:
        print("Не удалось загрузить данные!")
        return
    
    print(f"Загружено {len(df)} нулломеров")
    
    # Создание визуализаций
    create_comprehensive_visualizations(df, args.results_dir, args.k)
    create_network_visualization(df, args.results_dir)
    
    # Генерация отчета
    generate_summary_report(df, args.results_dir, args.k)
    
    print(f"\nВсе визуализации и отчеты сохранены в {args.results_dir}/")

if __name__ == '__main__':
    main()
