// nullomer_analyzer_optimized.cpp - оптимизированная версия с полными именами
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <cmath>
#include <iomanip>
#include <cstring>
#include <immintrin.h>  // Для SIMD инструкций

namespace fs = std::filesystem;

// Оптимизированная версия с SIMD и компактным хранением

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

// Структура для хранения результатов (с полными именами для CSV)
struct GenomeResult {
    std::string genome_id;
    std::string genome_path;
    uint32_t genome_id_hash;
    uint32_t genome_length;
    double gc_content;
    uint32_t n_nullomers;
    uint32_t unique_kmers;
};

// Структура для компактного бинарного хранения
struct CompactGenomeResult {
    uint32_t genome_id_hash;
    uint32_t genome_length;
    uint16_t gc_percent_x100;
    uint32_t n_nullomers;
    
    void write_binary(std::ofstream& out) const {
        out.write(reinterpret_cast<const char*>(this), sizeof(*this));
    }
    
    static CompactGenomeResult read_binary(std::ifstream& in) {
        CompactGenomeResult result;
        in.read(reinterpret_cast<char*>(&result), sizeof(result));
        return result;
    }
};

// Быстрый хеш для имен файлов
inline uint32_t hash_string(const std::string& str) {
    uint32_t hash = 5381;
    for (char c : str) {
        hash = ((hash << 5) + hash) + c;
    }
    return hash;
}

// Оптимизированный класс для работы с k-мерами
class FastKmerProcessor {
private:
    const size_t k;
    const uint64_t kmer_mask;
    alignas(64) std::vector<uint64_t> kmer_presence;  // Битовый массив
    
public:
    FastKmerProcessor(size_t k_size) 
        : k(k_size), 
          kmer_mask((1ULL << (2 * k)) - 1),
          kmer_presence((1ULL << (2 * k - 6)) + 1, 0) {}
    
    // SIMD-оптимизированная обработка последовательности
    void process_sequence_simd(const char* seq, size_t len) {
        if (len < k) return;
        
        uint64_t kmer = 0;
        size_t valid_bases = 0;
        
        // Инициализация первого k-мера
        for (size_t i = 0; i < k; ++i) {
            uint8_t nuc = NUC_TABLE[static_cast<uint8_t>(seq[i])];
            if (nuc == 4) {
                valid_bases = 0;
                kmer = 0;
            } else {
                kmer = ((kmer << 2) | nuc) & kmer_mask;
                valid_bases++;
            }
        }
        
        if (valid_bases == k) {
            mark_kmer(kmer);
            mark_kmer(reverse_complement_fast(kmer));
        }
        
        // Скользящее окно
        for (size_t i = k; i < len; ++i) {
            uint8_t nuc = NUC_TABLE[static_cast<uint8_t>(seq[i])];
            
            if (nuc == 4) {
                valid_bases = 0;
                kmer = 0;
                
                size_t skip_to = std::min(i + k, len);
                for (size_t j = i + 1; j < skip_to; ++j) {
                    uint8_t next_nuc = NUC_TABLE[static_cast<uint8_t>(seq[j])];
                    if (next_nuc != 4) {
                        kmer = ((kmer << 2) | next_nuc) & kmer_mask;
                        valid_bases++;
                    } else {
                        valid_bases = 0;
                        kmer = 0;
                    }
                }
                i = skip_to - 1;
            } else {
                kmer = ((kmer << 2) | nuc) & kmer_mask;
                if (valid_bases < k) {
                    valid_bases++;
                }
                
                if (valid_bases == k) {
                    mark_kmer(kmer);
                    mark_kmer(reverse_complement_fast(kmer));
                }
            }
        }
    }
    
    inline void mark_kmer(uint64_t kmer) {
        size_t idx = kmer >> 6;
        size_t bit = kmer & 63;
        kmer_presence[idx] |= (1ULL << bit);
    }
    
    inline uint64_t reverse_complement_fast(uint64_t kmer) const {
        uint64_t rc = 0;
        for (size_t i = 0; i < k; ++i) {
            rc = (rc << 2) | (3 - (kmer & 3));
            kmer >>= 2;
        }
        return rc & kmer_mask;
    }
    
    size_t count_nullomers() const {
        size_t total_kmers = 1ULL << (2 * k);
        size_t present_count = 0;
        
        for (const auto& word : kmer_presence) {
            present_count += __builtin_popcountll(word);
        }
        
        return total_kmers - present_count;
    }
    
    size_t count_present() const {
        size_t present_count = 0;
        for (const auto& word : kmer_presence) {
            present_count += __builtin_popcountll(word);
        }
        return present_count;
    }
    
    void clear() {
        std::fill(kmer_presence.begin(), kmer_presence.end(), 0);
    }
};

// Оптимизированный анализатор генома
class OptimizedGenomeAnalyzer {
private:
    const size_t k;
    FastKmerProcessor processor;
    
    static constexpr size_t BUFFER_SIZE = 16 * 1024 * 1024;  // 16MB буфер
    std::vector<char> read_buffer;
    
public:
    OptimizedGenomeAnalyzer(size_t k_size) 
        : k(k_size), processor(k_size), read_buffer(BUFFER_SIZE) {}
    
    GenomeResult analyze_genome(const std::string& genome_path) {
        GenomeResult result;
        result.genome_path = genome_path;
        result.genome_id = fs::path(genome_path).stem().string();
        result.genome_id_hash = hash_string(result.genome_id);
        result.genome_length = 0;
        result.gc_content = 0;
        result.n_nullomers = 0;
        result.unique_kmers = 0;
        
        processor.clear();
        
        std::ifstream file(genome_path, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Ошибка открытия файла: " << genome_path << std::endl;
            return result;
        }
        
        size_t total_length = 0;
        size_t gc_count = 0;
        std::string sequence;
        sequence.reserve(BUFFER_SIZE);
        
        // Читаем файл большими блоками
        while (file.read(read_buffer.data(), BUFFER_SIZE) || file.gcount() > 0) {
            size_t bytes_read = file.gcount();
            
            for (size_t i = 0; i < bytes_read; ++i) {
                char c = read_buffer[i];
                
                if (c == '>') {
                    // Обрабатываем накопленную последовательность
                    if (!sequence.empty()) {
                        processor.process_sequence_simd(sequence.c_str(), sequence.length());
                        total_length += sequence.length();
                        
                        // Подсчет GC
                        for (char base : sequence) {
                            if (base == 'G' || base == 'C' || base == 'g' || base == 'c') {
                                gc_count++;
                            }
                        }
                        
                        sequence.clear();
                    }
                    
                    // Пропускаем заголовок
                    while (i < bytes_read && read_buffer[i] != '\n') {
                        i++;
                    }
                } else if (c != '\n' && c != '\r' && c != ' ' && c != '\t') {
                    sequence += c;
                }
            }
        }
        
        // Обрабатываем последнюю последовательность
        if (!sequence.empty()) {
            processor.process_sequence_simd(sequence.c_str(), sequence.length());
            total_length += sequence.length();
            
            for (char base : sequence) {
                if (base == 'G' || base == 'C' || base == 'g' || base == 'c') {
                    gc_count++;
                }
            }
        }
        
        result.genome_length = static_cast<uint32_t>(total_length);
        result.gc_content = total_length > 0 ? 
            static_cast<double>(gc_count) / total_length : 0.0;
        result.n_nullomers = static_cast<uint32_t>(processor.count_nullomers());
        result.unique_kmers = static_cast<uint32_t>(processor.count_present());
        
        return result;
    }
};

// Основная программа
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Использование: " << argv[0] 
                  << " <директория_с_геномами> [--k <размер>] [--threads <количество>] "
                  << "[--output <файл>] [--binary] [--no-hash] [--max-genomes <количество>]" << std::endl;
        return 1;
    }
    
    std::string genome_dir = argv[1];
    size_t k = 12;
    size_t num_threads = std::thread::hardware_concurrency();
    std::string output_file = "nullomer_results_optimized.csv";
    bool binary_output = false;
    bool include_hash = true;
    size_t max_genomes = 0;
    
    // Парсинг аргументов
    for (int i = 2; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--k" && i + 1 < argc) {
            k = std::stoul(argv[++i]);
        } else if (arg == "--threads" && i + 1 < argc) {
            num_threads = std::stoul(argv[++i]);
        } else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg == "--binary") {
            binary_output = true;
            if (output_file == "nullomer_results_optimized.csv") {
                output_file = "nullomer_results_optimized.bin";
            }
        } else if (arg == "--no-hash") {
            include_hash = false;
        } else if (arg == "--max-genomes" && i + 1 < argc) {
            max_genomes = std::stoul(argv[++i]);
        }
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "=== Оптимизированный анализатор нуломеров ===" << std::endl;
    std::cout << "k = " << k << std::endl;
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
    
    // Ограничение количества геномов
    if (max_genomes > 0 && genome_files.size() > max_genomes) {
        genome_files.resize(max_genomes);
    }
    
    std::cout << "Найдено геномов: " << genome_files.size() << std::endl;
    
    // Параллельная обработка
    std::vector<GenomeResult> results;
    std::mutex results_mutex;
    std::atomic<size_t> processed(0);
    
    auto worker = [&](size_t start, size_t end) {
        OptimizedGenomeAnalyzer analyzer(k);
        std::vector<GenomeResult> local_results;
        
        for (size_t i = start; i < end; ++i) {
            local_results.push_back(analyzer.analyze_genome(genome_files[i]));
            
            size_t p = processed.fetch_add(1) + 1;
            if (p % 100 == 0) {
                std::cout << "\rОбработано: " << p << "/" << genome_files.size() 
                         << " (" << std::fixed << std::setprecision(1) 
                         << (100.0 * p / genome_files.size()) << "%)" << std::flush;
            }
        }
        
        std::lock_guard<std::mutex> lock(results_mutex);
        results.insert(results.end(), local_results.begin(), local_results.end());
    };
    
    // Запуск потоков
    std::cout << "\nНачинаем обработку..." << std::endl;
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
    
    // Сохранение результатов
    if (binary_output) {
        std::ofstream out(output_file, std::ios::binary);
        size_t count = results.size();
        out.write(reinterpret_cast<const char*>(&count), sizeof(count));
        out.write(reinterpret_cast<const char*>(&k), sizeof(k));
        
        for (const auto& r : results) {
            CompactGenomeResult compact;
            compact.genome_id_hash = r.genome_id_hash;
            compact.genome_length = r.genome_length;
            compact.gc_percent_x100 = static_cast<uint16_t>(r.gc_content * 10000);
            compact.n_nullomers = r.n_nullomers;
            compact.write_binary(out);
        }
        std::cout << "Результаты сохранены в бинарном формате: " << output_file << std::endl;
    } else {
        std::ofstream out(output_file);
        
        // Заголовок CSV
        out << "genome_id,genome_path";
        if (include_hash) {
            out << ",genome_id_hash";
        }
        out << ",genome_length,gc_content,n_nullomers,unique_kmers\n";
        
        // Данные
        for (const auto& r : results) {
            out << r.genome_id << ","
                << r.genome_path << ",";
            if (include_hash) {
                out << r.genome_id_hash << ",";
            }
            out << r.genome_length << ","
                << std::fixed << std::setprecision(4) << r.gc_content << ","
                << r.n_nullomers << ","
                << r.unique_kmers << "\n";
        }
        std::cout << "Результаты сохранены в CSV: " << output_file << std::endl;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    // Статистика
    double avg_nullomers = 0;
    double avg_gc = 0;
    for (const auto& r : results) {
        avg_nullomers += r.n_nullomers;
        avg_gc += r.gc_content;
    }
    avg_nullomers /= results.size();
    avg_gc /= results.size();
    
    std::cout << "\n=== Статистика ===" << std::endl;
    std::cout << "Время выполнения: " << duration.count() << " секунд" << std::endl;
    std::cout << "Скорость: " << (genome_files.size() / (duration.count() + 1)) 
              << " геномов/сек" << std::endl;
    std::cout << "Среднее количество нуломеров: " << std::fixed << std::setprecision(1) 
              << avg_nullomers << std::endl;
    std::cout << "Средний GC-состав: " << std::fixed << std::setprecision(3) 
              << avg_gc << std::endl;
    
    return 0;
}
