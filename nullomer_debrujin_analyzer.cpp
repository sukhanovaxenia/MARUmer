// nullomer_debruijn_analyzer.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <thread>
#include <mutex>
#include <atomic>
#include <filesystem>
#include <iomanip>
#include <cstring>
#include <immintrin.h>
#include <bitset>
#include <queue>

namespace fs = std::filesystem;

// Константы
constexpr size_t K = 13;
constexpr size_t TOTAL_KMERS = 1ULL << (2 * K);  // 4^13
constexpr size_t K_MINUS_1 = K - 1;

// Таблица для быстрого преобразования нуклеотидов
alignas(64) static constexpr uint8_t NUC_TABLE[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // A, C, G
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // T
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // a, c, g
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // t
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

static const char NUCS[] = "ACGT";

// Структура для хранения информации о k-мере
struct KmerInfo {
    uint64_t kmer;
    uint32_t count;
    bool is_present;
    
    KmerInfo() : kmer(0), count(0), is_present(false) {}
};

// Структура для анализа контекста нулломера
struct NullomerContext {
    std::string nullomer;
    uint64_t nullomer_encoded;
    std::vector<std::string> possible_predecessors;
    std::vector<std::string> possible_successors;
    uint32_t predecessor_count;
    uint32_t successor_count;
    double gc_content;
    std::string classification;
};

// Класс для работы с графом де Брюйна
class DeBruijnGraph {
private:
    // Битовый массив для быстрой проверки присутствия k-мера
    std::vector<uint64_t> kmer_presence;
    // Счетчики k-меров (только для присутствующих)
    std::unordered_map<uint64_t, uint32_t> kmer_counts;
    
    // Кодирование k-мера
    inline uint64_t encode_kmer(const char* seq, size_t len) const {
        uint64_t kmer = 0;
        for (size_t i = 0; i < len; ++i) {
            uint8_t nuc = NUC_TABLE[static_cast<uint8_t>(seq[i])];
            if (nuc == 4) return UINT64_MAX;
            kmer = (kmer << 2) | nuc;
        }
        return kmer;
    }
    
    // Декодирование k-мера
    std::string decode_kmer(uint64_t encoded, size_t len) const {
        std::string kmer(len, 'N');
        for (int i = len - 1; i >= 0; --i) {
            kmer[i] = NUCS[encoded & 3];
            encoded >>= 2;
        }
        return kmer;
    }
    
    // Получение canonical k-мера
    uint64_t canonical_kmer(uint64_t kmer) const {
        uint64_t rc = reverse_complement(kmer);
        return std::min(kmer, rc);
    }
    
    // Reverse complement
    uint64_t reverse_complement(uint64_t kmer) const {
        uint64_t rc = 0;
        for (size_t i = 0; i < K; ++i) {
            rc = (rc << 2) | (3 - (kmer & 3));
            kmer >>= 2;
        }
        return rc & ((1ULL << (2 * K)) - 1);
    }
    
public:
    DeBruijnGraph() : kmer_presence((TOTAL_KMERS + 63) / 64, 0) {}
    
    // Пометить k-мер как присутствующий
    void mark_kmer(uint64_t kmer) {
        if (kmer >= TOTAL_KMERS) return;
        
        uint64_t canonical = canonical_kmer(kmer);
        size_t idx = canonical >> 6;
        size_t bit = canonical & 63;
        
        // Атомарная операция для потокобезопасности
        __sync_fetch_and_or(&kmer_presence[idx], 1ULL << bit);
        
        // Увеличиваем счетчик
        std::lock_guard<std::mutex> lock(count_mutex);
        kmer_counts[canonical]++;
    }
    
    // Проверка присутствия k-мера
    bool has_kmer(uint64_t kmer) const {
        if (kmer >= TOTAL_KMERS) return false;
        uint64_t canonical = canonical_kmer(kmer);
        size_t idx = canonical >> 6;
        size_t bit = canonical & 63;
        return kmer_presence[idx] & (1ULL << bit);
    }
    
    // Получение всех соседей k-мера в графе
    std::vector<uint64_t> get_neighbors(uint64_t kmer) const {
        std::vector<uint64_t> neighbors;
        
        // Суффикс длины k-1
        uint64_t suffix = kmer & ((1ULL << (2 * K_MINUS_1)) - 1);
        
        // Проверяем все возможные исходящие ребра
        for (int nuc = 0; nuc < 4; ++nuc) {
            uint64_t next_kmer = (suffix << 2) | nuc;
            if (has_kmer(next_kmer)) {
                neighbors.push_back(next_kmer);
            }
        }
        
        return neighbors;
    }
    
    // Получение всех предшественников k-мера
    std::vector<uint64_t> get_predecessors(uint64_t kmer) const {
        std::vector<uint64_t> predecessors;
        
        // Префикс длины k-1
        uint64_t prefix = kmer >> 2;
        
        // Проверяем все возможные входящие ребра
        for (int nuc = 0; nuc < 4; ++nuc) {
            uint64_t prev_kmer = (uint64_t(nuc) << (2 * K_MINUS_1)) | prefix;
            if (has_kmer(prev_kmer)) {
                predecessors.push_back(prev_kmer);
            }
        }
        
        return predecessors;
    }
    
    // Анализ контекста нулломера
    NullomerContext analyze_nullomer_context(uint64_t nullomer) const {
        NullomerContext context;
        context.nullomer_encoded = nullomer;
        context.nullomer = decode_kmer(nullomer, K);
        
        // Получаем потенциальных предшественников
        auto preds = get_predecessors(nullomer);
        context.predecessor_count = preds.size();
        for (auto pred : preds) {
            context.possible_predecessors.push_back(decode_kmer(pred, K));
        }
        
        // Получаем потенциальных последователей
        auto succs = get_neighbors(nullomer);
        context.successor_count = succs.size();
        for (auto succ : succs) {
            context.possible_successors.push_back(decode_kmer(succ, K));
        }
        
        // Вычисляем GC-содержание
        int gc_count = 0;
        for (char c : context.nullomer) {
            if (c == 'G' || c == 'C') gc_count++;
        }
        context.gc_content = static_cast<double>(gc_count) / K;
        
        // Классификация на основе контекста
        if (context.predecessor_count == 0 && context.successor_count == 0) {
            context.classification = "isolated";
        } else if (context.predecessor_count == 0) {
            context.classification = "no_predecessors";
        } else if (context.successor_count == 0) {
            context.classification = "no_successors";
        } else if (context.predecessor_count == 1 && context.successor_count == 1) {
            context.classification = "linear";
        } else {
            context.classification = "branching";
        }
        
        return context;
    }
    
    // Получение всех нулломеров
    std::vector<uint64_t> get_nullomers() const {
        std::vector<uint64_t> nullomers;
        
        for (uint64_t kmer = 0; kmer < TOTAL_KMERS; ++kmer) {
            if (!has_kmer(kmer)) {
                // Проверяем canonical форму
                uint64_t canonical = canonical_kmer(kmer);
                if (kmer == canonical) {  // Добавляем только canonical формы
                    nullomers.push_back(kmer);
                }
            }
        }
        
        return nullomers;
    }
    
    // Статистика графа
    void print_statistics() const {
        size_t present_count = 0;
        for (const auto& word : kmer_presence) {
            present_count += __builtin_popcountll(word);
        }
        
        std::cout << "=== Статистика графа де Брюйна ===" << std::endl;
        std::cout << "Всего возможных k-меров: " << TOTAL_KMERS << std::endl;
        std::cout << "Присутствует k-меров: " << present_count << std::endl;
        std::cout << "Нулломеров: " << (TOTAL_KMERS - present_count) << std::endl;
        std::cout << "Покрытие: " << std::fixed << std::setprecision(2) 
                  << (100.0 * present_count / TOTAL_KMERS) << "%" << std::endl;
    }
    
private:
    mutable std::mutex count_mutex;
};

// Класс для обработки геномов
class GenomeProcessor {
private:
    DeBruijnGraph& graph;
    static constexpr size_t BUFFER_SIZE = 16 * 1024 * 1024;  // 16MB буфер
    
public:
    GenomeProcessor(DeBruijnGraph& g) : graph(g) {}
    
    void process_genome(const std::string& genome_path) {
        std::ifstream file(genome_path, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Ошибка открытия файла: " << genome_path << std::endl;
            return;
        }
        
        std::vector<char> buffer(BUFFER_SIZE);
        std::string sequence;
        sequence.reserve(BUFFER_SIZE);
        
        while (file.read(buffer.data(), BUFFER_SIZE) || file.gcount() > 0) {
            size_t bytes_read = file.gcount();
            
            for (size_t i = 0; i < bytes_read; ++i) {
                char c = buffer[i];
                
                if (c == '>') {
                    // Обрабатываем накопленную последовательность
                    if (!sequence.empty()) {
                        process_sequence(sequence);
                        sequence.clear();
                    }
                    
                    // Пропускаем заголовок
                    while (i < bytes_read && buffer[i] != '\n') {
                        i++;
                    }
                } else if (c != '\n' && c != '\r' && c != ' ' && c != '\t') {
                    sequence += c;
                }
            }
        }
        
        // Обрабатываем последнюю последовательность
        if (!sequence.empty()) {
            process_sequence(sequence);
        }
    }
    
private:
    void process_sequence(const std::string& seq) {
        if (seq.length() < K) return;
        
        uint64_t kmer = 0;
        size_t valid_bases = 0;
        
        // Инициализация первого k-мера
        for (size_t i = 0; i < K; ++i) {
            uint8_t nuc = NUC_TABLE[static_cast<uint8_t>(seq[i])];
            if (nuc == 4) {
                valid_bases = 0;
                kmer = 0;
            } else {
                kmer = ((kmer << 2) | nuc) & ((1ULL << (2 * K)) - 1);
                valid_bases++;
            }
        }
        
        if (valid_bases == K) {
            graph.mark_kmer(kmer);
        }
        
        // Скользящее окно
        for (size_t i = K; i < seq.length(); ++i) {
            uint8_t nuc = NUC_TABLE[static_cast<uint8_t>(seq[i])];
            
            if (nuc == 4) {
                valid_bases = 0;
                kmer = 0;
            } else {
                kmer = ((kmer << 2) | nuc) & ((1ULL << (2 * K)) - 1);
                if (valid_bases < K) {
                    valid_bases++;
                }
                
                if (valid_bases == K) {
                    graph.mark_kmer(kmer);
                }
            }
        }
    }
};

// Анализатор структуры графа
class GraphStructureAnalyzer {
private:
    const DeBruijnGraph& graph;
    
public:
    GraphStructureAnalyzer(const DeBruijnGraph& g) : graph(g) {}
    
    // Анализ связанных компонент
    void analyze_connected_components() {
        std::cout << "\n=== Анализ связанных компонент ===" << std::endl;
        
        std::unordered_set<uint64_t> visited;
        std::vector<size_t> component_sizes;
        
        for (uint64_t kmer = 0; kmer < TOTAL_KMERS; ++kmer) {
            if (graph.has_kmer(kmer) && visited.find(kmer) == visited.end()) {
                size_t component_size = bfs_component(kmer, visited);
                component_sizes.push_back(component_size);
            }
        }
        
        std::sort(component_sizes.rbegin(), component_sizes.rend());
        
        std::cout << "Количество компонент: " << component_sizes.size() << std::endl;
        if (!component_sizes.empty()) {
            std::cout << "Размер наибольшей компоненты: " << component_sizes[0] << std::endl;
            std::cout << "Размер наименьшей компоненты: " << component_sizes.back() << std::endl;
        }
    }
    
    // Анализ плотности графа
    void analyze_graph_density() {
        std::cout << "\n=== Анализ плотности графа ===" << std::endl;
        
        std::vector<int> in_degrees(5, 0);   // 0, 1, 2, 3, 4 входящих ребер
        std::vector<int> out_degrees(5, 0);  // 0, 1, 2, 3, 4 исходящих ребер
        
        size_t total_edges = 0;
        
        for (uint64_t kmer = 0; kmer < TOTAL_KMERS; ++kmer) {
            if (graph.has_kmer(kmer)) {
                auto preds = graph.get_predecessors(kmer);
                auto succs = graph.get_neighbors(kmer);
                
                in_degrees[std::min(preds.size(), size_t(4))]++;
                out_degrees[std::min(succs.size(), size_t(4))]++;
                
                total_edges += succs.size();
            }
        }
        
        std::cout << "Распределение входящих степеней:" << std::endl;
        for (int i = 0; i <= 4; ++i) {
            std::cout << "  Степень " << i << ": " << in_degrees[i] << " узлов" << std::endl;
        }
        
        std::cout << "\nРаспределение исходящих степеней:" << std::endl;
        for (int i = 0; i <= 4; ++i) {
            std::cout << "  Степень " << i << ": " << out_degrees[i] << " узлов" << std::endl;
        }
        
        std::cout << "\nОбщее количество ребер: " << total_edges << std::endl;
    }
    
private:
    size_t bfs_component(uint64_t start, std::unordered_set<uint64_t>& visited) {
        std::queue<uint64_t> queue;
        queue.push(start);
        visited.insert(start);
        
        size_t component_size = 1;
        
        while (!queue.empty()) {
            uint64_t current = queue.front();
            queue.pop();
            
            // Проверяем всех соседей
            auto neighbors = graph.get_neighbors(current);
            for (auto neighbor : neighbors) {
                if (visited.find(neighbor) == visited.end()) {
                    visited.insert(neighbor);
                    queue.push(neighbor);
                    component_size++;
                }
            }
            
            // Также проверяем предшественников (для неориентированной связности)
            auto preds = graph.get_predecessors(current);
            for (auto pred : preds) {
                if (visited.find(pred) == visited.end()) {
                    visited.insert(pred);
                    queue.push(pred);
                    component_size++;
                }
            }
        }
        
        return component_size;
    }
};

// Основная функция
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Использование: " << argv[0] 
                  << " <директория_с_геномами> [--threads <количество>] "
                  << "[--output <файл>] [--max-genomes <количество>]" << std::endl;
        return 1;
    }
    
    std::string genome_dir = argv[1];
    size_t num_threads = std::thread::hardware_concurrency();
    std::string output_file = "nullomer_debruijn_analysis.txt";
    size_t max_genomes = 0;
    
    // Парсинг аргументов
    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--threads" && i + 1 < argc) {
            num_threads = std::stoul(argv[++i]);
        } else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg == "--max-genomes" && i + 1 < argc) {
            max_genomes = std::stoul(argv[++i]);
        }
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "=== Анализатор нулломеров с графом де Брюйна ===" << std::endl;
    std::cout << "k = " << K << std::endl;
    std::cout << "Потоков: " << num_threads << std::endl;
    std::cout << "Директория: " << genome_dir << std::endl;
    
    // Поиск файлов
    std::vector<std::string> genome_files;
    for (const auto& entry : fs::recursive_directory_iterator(genome_dir)) {
        if (entry.is_regular_file()) {
            std::string ext = entry.path().extension().string();
            std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
            if (ext == ".fna" || ext == ".fasta" || ext == ".fa") {
                genome_files.push_back(entry.path().string());
            }
        }
    }
    
    if (max_genomes > 0 && genome_files.size() > max_genomes) {
        genome_files.resize(max_genomes);
    }
    
    std::cout << "Найдено геномов: " << genome_files.size() << std::endl;
    
    // Создаем граф де Брюйна
    DeBruijnGraph graph;
    
    // Параллельная обработка геномов
    std::atomic<size_t> processed(0);
    std::mutex progress_mutex;
    
    auto worker = [&](size_t start, size_t end) {
        GenomeProcessor processor(graph);
        
        for (size_t i = start; i < end; ++i) {
            processor.process_genome(genome_files[i]);
            
            size_t p = processed.fetch_add(1) + 1;
            if (p % 100 == 0) {
                std::lock_guard<std::mutex> lock(progress_mutex);
                std::cout << "\rОбработано: " << p << "/" << genome_files.size() 
                         << " (" << std::fixed << std::setprecision(1) 
                         << (100.0 * p / genome_files.size()) << "%)" << std::flush;
            }
        }
    };
    
    // Запуск потоков
    std::cout << "\nПостроение графа де Брюйна..." << std::endl;
    std::vector<std::thread> threads;
    size_t chunk_size = (genome_files.size() + num_threads - 1) / num_threads;
    
    for (size_t i = 0; i < num_threads; ++i) {
        size_t start = i * chunk_size;
        size_t end = std::min(start + chunk_size, genome_files.size());
        if (start < end) {
            threads.emplace_back(worker, start, end);
        }
    }
    
    for (auto& t : threads) {
        t.join();
    }
    
    std::cout << "\rОбработано: " << genome_files.size() << "/" << genome_files.size() 
              << " (100.0%)" << std::endl;
    
    // Анализ графа
    graph.print_statistics();
    
    // Структурный анализ
    GraphStructureAnalyzer analyzer(graph);
    analyzer.analyze_graph_density();
    analyzer.analyze_connected_components();
    
    // Анализ нулломеров
    std::cout << "\n=== Анализ нулломеров ===" << std::endl;
    auto nullomers = graph.get_nullomers();
    std::cout << "Найдено нулломеров (canonical): " << nullomers.size() << std::endl;
    
    // Классификация нулломеров
    std::unordered_map<std::string, std::vector<NullomerContext>> classified_nullomers;
    
    for (auto nullomer : nullomers) {
        auto context = graph.analyze_nullomer_context(nullomer);
        classified_nullomers[context.classification].push_back(context);
    }
    
    // Сохранение результатов
    std::ofstream out(output_file);
    if (!out.is_open()) {
        std::cerr << "Ошибка создания выходного файла: " << output_file << std::endl;
        return 1;
    }
    
    out << "# Анализ нулломеров в контексте графа де Брюйна\n";
    out << "# k = " << K << "\n";
    out << "# Геномов обработано: " << genome_files.size() << "\n";
    out << "# Всего нулломеров: " << nullomers.size() << "\n\n";
    
    // Сводка по классификации
    out << "## Классификация нулломеров\n";
    for (const auto& [classification, contexts] : classified_nullomers) {
        out << classification << ": " << contexts.size() << " нулломеров\n";
    }
    out << "\n";
    
    // Детальная информация о нулломерах
    out << "## Детальный анализ нулломеров\n";
    out << "nullomer\tclassification\tgc_content\tpred_count\tsucc_count\tpredecessors\tsuccessors\n";
    
    for (const auto& [classification, contexts] : classified_nullomers) {
        for (const auto& ctx : contexts) {
            out << ctx.nullomer << "\t"
                << ctx.classification << "\t"
                << std::fixed << std::setprecision(3) << ctx.gc_content << "\t"
                << ctx.predecessor_count << "\t"
                << ctx.successor_count << "\t";
            
            // Предшественники
            if (ctx.possible_predecessors.empty()) {
                out << "none";
            } else {
                for (size_t i = 0; i < ctx.possible_predecessors.size(); ++i) {
                    if (i > 0) out << ",";
                    out << ctx.possible_predecessors[i];
                }
            }
            out << "\t";
            
            // Последователи
            if (ctx.possible_successors.empty()) {
                out << "none";
            } else {
                for (size_t i = 0; i < ctx.possible_successors.size(); ++i) {
                    if (i > 0) out << ",";
                    out << ctx.possible_successors[i];
                }
            }
            out << "\n";
        }
    }
    
    out.close();
    std::cout << "\nРезультаты сохранены в " << output_file << std::endl;
    
    // Статистика по GC-содержанию нулломеров
    std::cout << "\n=== Анализ GC-содержания нулломеров ===" << std::endl;
    std::vector<int> gc_distribution(K + 1, 0);
    
    for (auto nullomer : nullomers) {
        std::string seq = graph.decode_kmer(nullomer, K);
        int gc_count = 0;
        for (char c : seq) {
            if (c == 'G' || c == 'C') gc_count++;
        }
        gc_distribution[gc_count]++;
    }
    
    std::cout << "Распределение GC-содержания:" << std::endl;
    for (int i = 0; i <= K; ++i) {
        if (gc_distribution[i] > 0) {
            double gc_percent = 100.0 * i / K;
            std::cout << "  GC=" << std::fixed << std::setprecision(1) << gc_percent 
                     << "%: " << gc_distribution[i] << " нулломеров" << std::endl;
        }
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    std::cout << "\nВремя выполнения: " << duration.count() << " секунд" << std::endl;
    
    return 0;
}
