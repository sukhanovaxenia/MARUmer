// core_nullomer_finder.cpp
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <thread>
#include <mutex>
#include <atomic>
#include <filesystem>
#include <iomanip>

namespace fs = std::filesystem;

// Класс для эффективного подсчета core нуломеров
class CoreNullomerFinder {
private:
    const size_t k;
    const size_t total_kmers;
    const double threshold;
    std::vector<std::atomic<uint32_t>> nullomer_counts;  // Счетчики для каждого k-мера
    
    // Кодирование k-мера
    uint64_t encode_kmer(const std::string& kmer) const {
        uint64_t encoded = 0;
        for (char c : kmer) {
            uint8_t nuc = 4;
            switch (c) {
                case 'A': case 'a': nuc = 0; break;
                case 'C': case 'c': nuc = 1; break;
                case 'G': case 'g': nuc = 2; break;
                case 'T': case 't': nuc = 3; break;
            }
            if (nuc == 4) return UINT64_MAX;
            encoded = (encoded << 2) | nuc;
        }
        return encoded;
    }
    
    // Декодирование k-мера
    std::string decode_kmer(uint64_t encoded) const {
        static const char nucs[] = "ACGT";
        std::string kmer(k, 'N');
        for (int i = k - 1; i >= 0; --i) {
            kmer[i] = nucs[encoded & 3];
            encoded >>= 2;
        }
        return kmer;
    }
    
public:
    CoreNullomerFinder(size_t k_size, double core_threshold) 
        : k(k_size), 
          total_kmers(1ULL << (2 * k)),
          threshold(core_threshold),
          nullomer_counts(total_kmers) {
        // Инициализируем все счетчики нулями
        for (auto& count : nullomer_counts) {
            count.store(0);
        }
    }
    
    // Обработка результатов одного генома
    void process_genome_result(const std::string& nullomers_line) {
        if (nullomers_line.empty()) return;
        
        // Создаем битовый массив для отметки присутствующих k-меров
        std::vector<bool> present(total_kmers, false);
        
        // Парсим строку с нуломерами (предполагаем формат: kmer1,kmer2,kmer3,...)
        size_t pos = 0;
        while (pos < nullomers_line.length()) {
            size_t comma_pos = nullomers_line.find(',', pos);
            if (comma_pos == std::string::npos) {
                comma_pos = nullomers_line.length();
            }
            
            std::string kmer = nullomers_line.substr(pos, comma_pos - pos);
            if (!kmer.empty()) {
                uint64_t encoded = encode_kmer(kmer);
                if (encoded != UINT64_MAX && encoded < total_kmers) {
                    present[encoded] = true;
                }
            }
            
            pos = comma_pos + 1;
        }
        
        // Увеличиваем счетчики для отсутствующих k-меров
        for (size_t i = 0; i < total_kmers; ++i) {
            if (!present[i]) {
                nullomer_counts[i].fetch_add(1);
            }
        }
    }
    
    // Загрузка из бинарного файла с битовыми массивами
    void process_binary_file(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) return;
        
        // Читаем битовый массив
        std::vector<uint8_t> bit_array((total_kmers + 7) / 8);
        file.read(reinterpret_cast<char*>(bit_array.data()), bit_array.size());
        
        // Обновляем счетчики
        for (size_t i = 0; i < total_kmers; ++i) {
            bool is_present = bit_array[i / 8] & (1 << (i % 8));
            if (!is_present) {
                nullomer_counts[i].fetch_add(1);
            }
        }
    }
    
    // Получение core нуломеров
    std::vector<std::pair<std::string, double>> get_core_nullomers(size_t total_genomes) {
        std::vector<std::pair<std::string, double>> core_nullomers;
        size_t threshold_count = static_cast<size_t>(total_genomes * threshold);
        
        for (size_t i = 0; i < total_kmers; ++i) {
            uint32_t count = nullomer_counts[i].load();
            if (count >= threshold_count) {
                std::string kmer = decode_kmer(i);
                double frequency = static_cast<double>(count) / total_genomes;
                core_nullomers.emplace_back(kmer, frequency);
            }
        }
        
        // Сортируем по частоте
        std::sort(core_nullomers.begin(), core_nullomers.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });
        
        return core_nullomers;
    }
    
    // Параллельная обработка CSV файла
    void process_csv_file(const std::string& csv_file, size_t num_threads = 4) {
        std::ifstream file(csv_file);
        if (!file.is_open()) {
            std::cerr << "Не удалось открыть файл: " << csv_file << std::endl;
            return;
        }
        
        // Пропускаем заголовок
        std::string header;
        std::getline(file, header);
        
        // Читаем все строки
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(file, line)) {
            lines.push_back(line);
        }
        file.close();
        
        // Обрабатываем параллельно
        std::atomic<size_t> current_line(0);
        std::vector<std::thread> threads;
        
        auto worker = [&]() {
            while (true) {
                size_t idx = current_line.fetch_add(1);
                if (idx >= lines.size()) break;
                
                // Простой парсер CSV - извлекаем последнее поле (nullomers)
                const std::string& line = lines[idx];
                size_t last_comma = line.rfind(',');
                if (last_comma != std::string::npos) {
                    std::string nullomers = line.substr(last_comma + 1);
                    process_genome_result(nullomers);
                }
            }
        };
        
        for (size_t i = 0; i < num_threads; ++i) {
            threads.emplace_back(worker);
        }
        
        for (auto& t : threads) {
            t.join();
        }
    }
};

// Основная функция для поиска core нуломеров
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Использование: " << argv[0] 
                  << " <результаты.csv> [--k <размер>] [--threshold <порог>] "
                  << "[--output <файл>] [--threads <количество>]" << std::endl;
        return 1;
    }
    
    std::string input_file = argv[1];
    size_t k = 12;
    double threshold = 0.95;
    std::string output_file = "core_nullomers.txt";
    size_t num_threads = std::thread::hardware_concurrency();
    
    // Парсинг аргументов
    for (int i = 2; i < argc; i += 2) {
        if (i + 1 >= argc) break;
        
        std::string arg = argv[i];
        if (arg == "--k") {
            k = std::stoul(argv[i + 1]);
        } else if (arg == "--threshold") {
            threshold = std::stod(argv[i + 1]);
        } else if (arg == "--output") {
            output_file = argv[i + 1];
        } else if (arg == "--threads") {
            num_threads = std::stoul(argv[i + 1]);
        }
    }
    
    std::cout << "=== Поиск Core Нуломеров ===" << std::endl;
    std::cout << "k = " << k << std::endl;
    std::cout << "Порог = " << (threshold * 100) << "%" << std::endl;
    std::cout << "Потоков = " << num_threads << std::endl;
    
    // Создаем поисковик
    CoreNullomerFinder finder(k, threshold);
    
    // Определяем тип входного файла
    if (input_file.substr(input_file.find_last_of(".") + 1) == "csv") {
        std::cout << "Обработка CSV файла..." << std::endl;
        finder.process_csv_file(input_file, num_threads);
    } else {
        std::cerr << "Неподдерживаемый формат файла" << std::endl;
        return 1;
    }
    
    // Подсчитываем количество геномов из CSV
    std::ifstream count_file(input_file);
    size_t genome_count = 0;
    std::string line;
    while (std::getline(count_file, line)) {
        genome_count++;
    }
    genome_count--; // Вычитаем заголовок
    count_file.close();
    
    std::cout << "Обработано геномов: " << genome_count << std::endl;
    
    // Получаем core нуломеры
    std::cout << "Поиск core нуломеров..." << std::endl;
    auto core_nullomers = finder.get_core_nullomers(genome_count);
    
    std::cout << "Найдено core нуломеров: " << core_nullomers.size() << std::endl;
    
    // Сохраняем результаты
    std::ofstream out(output_file);
    if (!out.is_open()) {
        std::cerr << "Не удалось создать выходной файл: " << output_file << std::endl;
        return 1;
    }
    
    out << "# Core nullomers (отсутствуют в >= " << (threshold * 100) << "% геномов)" << std::endl;
    out << "# kmer\tfrequency" << std::endl;
    
    for (const auto& [kmer, freq] : core_nullomers) {
        out << kmer << "\t" << std::fixed << std::setprecision(3) << freq << std::endl;
    }
    
    out.close();
    std::cout << "Результаты сохранены в " << output_file << std::endl;
    
    // Выводим примеры
    if (!core_nullomers.empty()) {
        std::cout << "\nПримеры core нуломеров:" << std::endl;
        size_t examples = std::min(size_t(10), core_nullomers.size());
        for (size_t i = 0; i < examples; ++i) {
            std::cout << core_nullomers[i].first << "\t" 
                     << std::fixed << std::setprecision(3) 
                     << core_nullomers[i].second << std::endl;
        }
    }
    
    return 0;
}
