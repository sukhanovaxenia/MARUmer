# Makefile для пайплайна анализа нулломеров

CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -pthread -fopenmp -Wall
LDFLAGS = -pthread -fopenmp

# Исполняемые файлы
TARGETS = nullomer_analyzer nullomer_analyzer_optimized nullomer_debruijn_analyzer core_nullomer_finder

# Исходные файлы
SOURCES = nullomer_analyzer.cpp nullomer_analyzer_optimized.cpp nullomer_debruijn_analyzer.cpp core_nullomer_finder.cpp

# Директории
BUILD_DIR = build
BIN_DIR = bin

# Python зависимости
PYTHON_DEPS = numpy pandas matplotlib seaborn scipy scikit-learn networkx tqdm

.PHONY: all clean install-deps test run-example

all: directories $(TARGETS)

directories:
	@mkdir -p $(BUILD_DIR) $(BIN_DIR)

# Компиляция основных программ
nullomer_analyzer: nullomer_analyzer.cpp
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@ $(LDFLAGS)

nullomer_analyzer_optimized: nullomer_analyzer_optimized.cpp
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@ $(LDFLAGS)

nullomer_debruijn_analyzer: nullomer_debruijn_analyzer.cpp
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@ $(LDFLAGS)

core_nullomer_finder: core_nullomer_finder.cpp
	$(CXX) $(CXXFLAGS) $< -o $(BIN_DIR)/$@ $(LDFLAGS)

# Установка Python зависимостей
install-deps:
	@echo "Установка Python зависимостей..."
	pip install $(PYTHON_DEPS)
	@echo "Проверка наличия Jellyfish..."
	@which jellyfish > /dev/null || echo "ВНИМАНИЕ: Jellyfish не установлен. Установите с помощью: sudo apt-get install jellyfish"

# Очистка
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)
	rm -f nullomer_analyzer nullomer_analyzer_optimized nullomer_debruijn_analyzer core_nullomer_finder

# Тестовый запуск
test: all
	@echo "Запуск тестов..."
	@if [ -d "test_genomes" ]; then \
		./run_nullomer_pipeline.sh -d test_genomes -k 13 -t 4 --max-genomes 10; \
	else \
		echo "Директория test_genomes не найдена. Создайте её и добавьте тестовые геномы."; \
	fi

# Пример использования
run-example:
	@echo "=== Пример запуска пайплайна ==="
	@echo ""
	@echo "1. Базовый C++ анализ:"
	@echo "   ./run_nullomer_pipeline.sh -d /path/to/genomes -m cpp"
	@echo ""
	@echo "2. Анализ с Jellyfish:"
	@echo "   ./run_nullomer_pipeline.sh -d /path/to/genomes -m jellyfish"
	@echo ""
	@echo "3. Анализ подмножества геномов:"
	@echo "   ./run_nullomer_pipeline.sh -d /path/to/genomes --max-genomes 1000"
	@echo ""
	@echo "4. Прямой запуск анализатора де Брюйна:"
	@echo "   ./bin/nullomer_debruijn_analyzer /path/to/genomes --threads 16"

# Создание дистрибутива
dist: all
	@echo "Создание дистрибутива..."
	@mkdir -p nullomer_pipeline_dist
	@cp -r bin nullomer_pipeline_dist/
	@cp *.py *.sh Makefile README.md nullomer_pipeline_dist/
	@tar -czf nullomer_pipeline_v1.0.tar.gz nullomer_pipeline_dist
	@rm -rf nullomer_pipeline_dist
	@echo "Дистрибутив создан: nullomer_pipeline_v1.0.tar.gz"

# Документация
docs:
	@echo "Генерация документации..."
	@mkdir -p docs
	@echo "# Документация пайплайна анализа нулломеров" > docs/README.md
	@echo "" >> docs/README.md
	@echo "## Компоненты:" >> docs/README.md
	@echo "- nullomer_analyzer: Базовый анализатор" >> docs/README.md
	@echo "- nullomer_analyzer_optimized: Оптимизированная версия с SIMD" >> docs/README.md
	@echo "- nullomer_debruijn_analyzer: Анализ с графом де Брюйна" >> docs/README.md
	@echo "- core_nullomer_finder: Поиск core нулломеров" >> docs/README.md
