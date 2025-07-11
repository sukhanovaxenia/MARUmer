#!/usr/bin/env python3
"""
analyze_jellyfish_output.py - Анализ вывода Jellyfish для поиска нулломеров
и построения графа де Брюйна
"""

import argparse
import sys
from collections import defaultdict
from itertools import product
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import networkx as nx
import json

def reverse_complement(seq):
    """Получение обратно-комплементарной последовательности"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

def canonical_kmer(kmer):
    """Получение канонической формы k-мера"""
    rc = reverse_complement(kmer)
    return min(kmer, rc)

def get_all_kmers(k):
    """Генерация всех возможных k-меров"""
    return [''.join(bases) for bases in product('ACGT', repeat=k)]

def load_jellyfish_kmers(filename):
    """Загрузка k-меров из вывода Jellyfish"""
    observed_kmers = set()
    kmer_counts = {}
    
    print(f"Загрузка k-меров из {filename}...")
    with open(filename, 'r') as f:
        for line in tqdm(f, desc="Чтение k-меров"):
            line = line.strip()
            if line.startswith('>'):
                count = int(line[1:])
            else:
                kmer = line.upper()
                canonical = canonical_kmer(kmer)
                observed_kmers.add(canonical)
                kmer_counts[canonical] = kmer_counts.get(canonical, 0) + count
    
    return observed_kmers, kmer_counts

def build_debruijn_graph(observed_kmers, k):
    """Построение графа де Брюйна из наблюдаемых k-меров"""
    G = nx.DiGraph()
    
    print("Построение графа де Брюйна...")
    for kmer in tqdm(observed_kmers, desc="Добавление узлов"):
        G.add_node(kmer)
        
        # Добавление ребер
        suffix = kmer[1:]
        for base in 'ACGT':
            next_kmer = suffix + base
            if canonical_kmer(next_kmer) in observed_kmers:
                G.add_edge(kmer, canonical_kmer(next_kmer))
    
    return G

def analyze_nullomer_context(nullomer, observed_kmers, k):
    """Анализ контекста нулломера в графе"""
    context = {
        'nullomer': nullomer,
        'canonical': canonical_kmer(nullomer),
        'predecessors': [],
        'successors': [],
        'gc_content': 0
    }
    
    # GC-содержание
    gc_count = sum(1 for base in nullomer if base in 'GC')
    context['gc_content'] = gc_count / k
    
    # Проверка возможных предшественников
    prefix = nullomer[:-1]
    for base in 'ACGT':
        pred_kmer = base + prefix
        if canonical_kmer(pred_kmer) in observed_kmers:
            context['predecessors'].append(pred_kmer)
    
    # Проверка возможных последователей
    suffix = nullomer[1:]
    for base in 'ACGT':
        succ_kmer = suffix + base
        if canonical_kmer(succ_kmer) in observed_kmers:
            context['successors'].append(succ_kmer)
    
    # Классификация
    if not context['predecessors'] and not context['successors']:
        context['classification'] = 'isolated'
    elif not context['predecessors']:
        context['classification'] = 'no_predecessors'
    elif not context['successors']:
        context['classification'] = 'no_successors'
    elif len(context['predecessors']) == 1 and len(context['successors']) == 1:
        context['classification'] = 'linear'
    else:
        context['classification'] = 'branching'
    
    return context

def analyze_graph_structure(G):
    """Анализ структуры графа де Брюйна"""
    print("\nАнализ структуры графа...")
    
    stats = {
        'nodes': G.number_of_nodes(),
        'edges': G.number_of_edges(),
        'density': nx.density(G),
        'components': nx.number_weakly_connected_components(G),
        'largest_component': len(max(nx.weakly_connected_components(G), key=len)) if G.number_of_nodes() > 0 else 0
    }
    
    # Распределение степеней
    in_degrees = dict(G.in_degree())
    out_degrees = dict(G.out_degree())
    
    stats['in_degree_dist'] = np.bincount([d for d in in_degrees.values()])
    stats['out_degree_dist'] = np.bincount([d for d in out_degrees.values()])
    
    return stats

def visualize_results(nullomer_contexts, graph_stats, output_dir, k):
    """Создание визуализаций результатов"""
    print("\nСоздание визуализаций...")
    
    # 1. Распределение классификаций нулломеров
    classifications = [ctx['classification'] for ctx in nullomer_contexts]
    class_counts = defaultdict(int)
    for c in classifications:
        class_counts[c] += 1
    
    plt.figure(figsize=(10, 6))
    plt.bar(class_counts.keys(), class_counts.values())
    plt.xlabel('Классификация нулломеров')
    plt.ylabel('Количество')
    plt.title(f'Распределение нулломеров по типам контекста (k={k})')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/nullomer_classifications.png')
    plt.close()
    
    # 2. GC-содержание нулломеров
    gc_contents = [ctx['gc_content'] for ctx in nullomer_contexts]
    
    plt.figure(figsize=(10, 6))
    plt.hist(gc_contents, bins=k+1, edgecolor='black', alpha=0.7)
    plt.xlabel('GC-содержание')
    plt.ylabel('Количество нулломеров')
    plt.title(f'Распределение GC-содержания нулломеров (k={k})')
    plt.axvline(x=0.5, color='red', linestyle='--', label='50% GC')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{output_dir}/nullomer_gc_distribution.png')
    plt.close()
    
    # 3. Распределение степеней в графе
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Входящие степени
    in_degrees = graph_stats['in_degree_dist']
    ax1.bar(range(len(in_degrees)), in_degrees)
    ax1.set_xlabel('Входящая степень')
    ax1.set_ylabel('Количество узлов')
    ax1.set_title('Распределение входящих степеней')
    
    # Исходящие степени
    out_degrees = graph_stats['out_degree_dist']
    ax2.bar(range(len(out_degrees)), out_degrees)
    ax2.set_xlabel('Исходящая степень')
    ax2.set_ylabel('Количество узлов')
    ax2.set_title('Распределение исходящих степеней')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/graph_degree_distribution.png')
    plt.close()
    
    # 4. Тепловая карта связности нулломеров
    connectivity_matrix = np.zeros((5, 5))  # max 4 predecessors/successors
    for ctx in nullomer_contexts:
        pred_count = min(len(ctx['predecessors']), 4)
        succ_count = min(len(ctx['successors']), 4)
        connectivity_matrix[pred_count, succ_count] += 1
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(connectivity_matrix, annot=True, fmt='g', cmap='YlOrRd',
                xticklabels=range(5), yticklabels=range(5))
    plt.xlabel('Количество последователей')
    plt.ylabel('Количество предшественников')
    plt.title('Матрица связности нулломеров')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/nullomer_connectivity_heatmap.png')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Анализ вывода Jellyfish для поиска нулломеров')
    parser.add_argument('jellyfish_output', help='Файл с выводом Jellyfish')
    parser.add_argument('--k', type=int, default=13, help='Размер k-мера')
    parser.add_argument('--output', default='nullomer_analysis.txt', help='Выходной файл')
    parser.add_argument('--visualize', action='store_true', help='Создать визуализации')
    parser.add_argument('--output-dir', default='.', help='Директория для сохранения результатов')
    
    args = parser.parse_args()
    
    # Загрузка наблюдаемых k-меров
    observed_kmers, kmer_counts = load_jellyfish_kmers(args.jellyfish_output)
    print(f"Загружено {len(observed_kmers)} уникальных k-меров")
    
    # Поиск нулломеров
    print("\nПоиск нулломеров...")
    all_kmers = set(canonical_kmer(kmer) for kmer in get_all_kmers(args.k))
    nullomers = all_kmers - observed_kmers
    print(f"Найдено {len(nullomers)} нулломеров")
    
    # Построение графа де Брюйна
    G = build_debruijn_graph(observed_kmers, args.k)
    graph_stats = analyze_graph_structure(G)
    
    # Анализ контекста нулломеров
    print("\nАнализ контекста нулломеров...")
    nullomer_contexts = []
    for nullomer in tqdm(nullomers, desc="Анализ нулломеров"):
        context = analyze_nullomer_context(nullomer, observed_kmers, args.k)
        nullomer_contexts.append(context)
    
    # Сортировка по классификации
    nullomer_contexts.sort(key=lambda x: (x['classification'], x['gc_content']))
    
    # Сохранение результатов
    with open(args.output, 'w') as f:
        # Заголовок
        f.write(f"# Анализ нулломеров (k={args.k})\n")
        f.write(f"# Всего k-меров: {len(all_kmers)}\n")
        f.write(f"# Наблюдаемых k-меров: {len(observed_kmers)}\n")
        f.write(f"# Нулломеров: {len(nullomers)}\n")
        f.write(f"# Покрытие: {100.0 * len(observed_kmers) / len(all_kmers):.2f}%\n\n")
        
        # Статистика графа
        f.write("## Статистика графа де Брюйна\n")
        f.write(f"Узлов: {graph_stats['nodes']}\n")
        f.write(f"Ребер: {graph_stats['edges']}\n")
        f.write(f"Плотность: {graph_stats['density']:.4f}\n")
        f.write(f"Компонент связности: {graph_stats['components']}\n")
        f.write(f"Размер наибольшей компоненты: {graph_stats['largest_component']}\n\n")
        
        # Сводка по классификациям
        f.write("## Классификация нулломеров\n")
        class_counts = defaultdict(int)
        for ctx in nullomer_contexts:
            class_counts[ctx['classification']] += 1
        
        for classification, count in sorted(class_counts.items()):
            f.write(f"{classification}: {count} ({100.0 * count / len(nullomers):.1f}%)\n")
        f.write("\n")
        
        # Детальная информация
        f.write("## Детальный список нулломеров\n")
        f.write("nullomer\tclassification\tgc_content\tpred_count\tsucc_count\tpredecessors\tsuccessors\n")
        
        for ctx in nullomer_contexts:
            predecessors = ','.join(ctx['predecessors']) if ctx['predecessors'] else 'none'
            successors = ','.join(ctx['successors']) if ctx['successors'] else 'none'
            
            f.write(f"{ctx['nullomer']}\t{ctx['classification']}\t{ctx['gc_content']:.3f}\t"
                   f"{len(ctx['predecessors'])}\t{len(ctx['successors'])}\t"
                   f"{predecessors}\t{successors}\n")
    
    print(f"\nРезультаты сохранены в {args.output}")
    
    # Сохранение JSON для дальнейшего анализа
    json_output = args.output.replace('.txt', '.json')
    with open(json_output, 'w') as f:
        json.dump({
            'k': args.k,
            'total_kmers': len(all_kmers),
            'observed_kmers': len(observed_kmers),
            'nullomers': len(nullomers),
            'graph_stats': {k: v if not isinstance(v, np.ndarray) else v.tolist() 
                           for k, v in graph_stats.items()},
            'nullomer_contexts': nullomer_contexts[:1000]  # Сохраняем первые 1000 для примера
        }, f, indent=2)
    
    # Создание визуализаций
    if args.visualize:
        visualize_results(nullomer_contexts, graph_stats, args.output_dir, args.k)
        print(f"Визуализации сохранены в {args.output_dir}/")

if __name__ == '__main__':
    main()
