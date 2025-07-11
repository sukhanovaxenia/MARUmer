#!/bin/bash
# run_nullomer_pipeline.sh - Полный пайплайн для анализа нулломеров с графом де Брюйна

# Цвета для вывода
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Параметры по умолчанию
GENOME_DIR=""
K=13
THREADS=$(nproc)
METHOD="cpp"  # cpp, jellyfish, aindex
OUTPUT_DIR="nullomer_analysis_results"
MAX_GENOMES=0
SKIP_BUILD=false

# Функция показа помощи
show_help() {
    echo -e "${GREEN}Пайплайн анализа нулломеров с графом де Брюйна${NC}"
    echo ""
    echo -e "${BLUE}Использование:${NC} $0 -d <genome_directory> [опции]"
    echo ""
    echo -e "${BLUE}Обязательные параметры:${NC}"
    echo "  -d, --dir PATH         Директория с геномами"
    echo ""
    echo -e "${BLUE}Опции:${NC}"
    echo "  -m, --method METHOD    Метод анализа: cpp (default), jellyfish, aindex"
    echo "  -k, --kmer SIZE        Размер k-мера (по умолчанию: 13)"
    echo "  -t, --threads N        Количество потоков (по умолчанию: $(nproc))"
    echo "  -o, --output DIR       Выходная директория (по умолчанию: nullomer_analysis_results)"
    echo "  --max-genomes N        Максимальное количество геномов для анализа"
    echo "  --skip-build           Пропустить компиляцию программ"
    echo "  -h, --help             Показать эту справку"
}

# Парсинг аргументов
while [[ $# -gt 0 ]]; do
    case $1 in
        -d|--dir)
            GENOME_DIR="$2"
            shift 2
            ;;
        -m|--method)
            METHOD="$2"
            shift 2
            ;;
        -k|--kmer)
            K="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --max-genomes)
            MAX_GENOMES="$2"
            shift 2
            ;;
        --skip-build)
            SKIP_BUILD=true
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo -e "${RED}Неизвестный аргумент: $1${NC}"
            show_help
            exit 1
            ;;
    esac
done

# Проверка обязательных параметров
if [ -z "$GENOME_DIR" ]; then
    echo -e "${RED}Ошибка: не указана директория с геномами${NC}"
    show_help
    exit 1
fi

if [ ! -d "$GENOME_DIR" ]; then
    echo -e "${RED}Ошибка: директория '$GENOME_DIR' не существует${NC}"
    exit 1
fi

# Создание выходной директории
mkdir -p "$OUTPUT_DIR"

# Функция для компиляции программ
compile_programs() {
    echo -e "${YELLOW}Компиляция программ...${NC}"
    
    # Компиляция анализатора с графом де Брюйна
    if [ ! -f "nullomer_debruijn_analyzer" ] || [ "$SKIP_BUILD" = false ]; then
        echo "Компиляция nullomer_debruijn_analyzer..."
        g++ -std=c++17 -O3 -march=native -pthread -fopenmp \
            nullomer_debruijn_analyzer.cpp -o nullomer_debruijn_analyzer
        
        if [ $? -ne 0 ]; then
            echo -e "${RED}Ошибка компиляции nullomer_debruijn_analyzer${NC}"
            exit 1
        fi
    fi
    
    echo -e "${GREEN}✓ Компиляция завершена${NC}"
}

# Функция для метода C++
run_cpp_method() {
    echo -e "${BLUE}Запуск C++ анализатора...${NC}"
    
    ./nullomer_debruijn_analyzer "$GENOME_DIR" \
        --threads $THREADS \
        --output "$OUTPUT_DIR/debruijn_analysis.txt" \
        $([ $MAX_GENOMES -gt 0 ] && echo "--max-genomes $MAX_GENOMES")
    
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ C++ анализ завершен${NC}"
    else
        echo -e "${RED}✗ Ошибка в C++ анализе${NC}"
        exit 1
    fi
}

# Функция для метода Jellyfish
run_jellyfish_method() {
    echo -e "${BLUE}Запуск Jellyfish метода...${NC}"
    
    # Проверка наличия Jellyfish
    if ! command -v jellyfish &> /dev/null; then
        echo -e "${RED}Ошибка: Jellyfish не установлен${NC}"
        echo "Установите Jellyfish: sudo apt-get install jellyfish"
        exit 1
    fi
    
    # Создание mega.fasta
    echo "Создание mega.fasta..."
    MEGAFASTA="$OUTPUT_DIR/mega.fasta"
    > "$MEGAFASTA"
    
    GENOME_COUNT=0
    for genome in "$GENOME_DIR"/*.{fna,fasta,fa}; do
        if [ -f "$genome" ]; then
            cat "$genome" >> "$MEGAFASTA"
            echo >> "$MEGAFASTA"  # Добавляем перевод строки между файлами
            
            GENOME_COUNT=$((GENOME_COUNT + 1))
            if [ $MAX_GENOMES -gt 0 ] && [ $GENOME_COUNT -ge $MAX_GENOMES ]; then
                break
            fi
        fi
    done
    
    echo "Объединено геномов: $GENOME_COUNT"
    
    # Запуск Jellyfish
    echo "Подсчет k-меров с Jellyfish..."
    jellyfish count -m $K -s 10G -t $THREADS -C "$MEGAFASTA" -o "$OUTPUT_DIR/kmers.jf"
    
    if [ $? -ne 0 ]; then
        echo -e "${RED}Ошибка в Jellyfish count${NC}"
        exit 1
    fi
    
    # Экспорт k-меров
    echo "Экспорт k-меров..."
    jellyfish dump "$OUTPUT_DIR/kmers.jf" > "$OUTPUT_DIR/kmers_${K}.txt"
    
    # Анализ с помощью Python скрипта
    echo "Анализ нулломеров..."
    python3 analyze_jellyfish_output.py \
        "$OUTPUT_DIR/kmers_${K}.txt" \
        --k $K \
        --output "$OUTPUT_DIR/jellyfish_nullomers.txt"
    
    echo -e "${GREEN}✓ Jellyfish анализ завершен${NC}"
}

# Функция для метода AIndex
run_aindex_method() {
    echo -e "${BLUE}Запуск AIndex метода...${NC}"
    
    # Проверка наличия AIndex
    if [ ! -f "./aindex" ]; then
        echo -e "${YELLOW}AIndex не найден, пытаемся скомпилировать...${NC}"
        if [ -d "Tools/aindex" ]; then
            cd Tools/aindex
            make
            cp aindex ../../
            cd ../..
        else
            echo -e "${RED}Ошибка: исходный код AIndex не найден${NC}"
            exit 1
        fi
    fi
    
    # Создание списка файлов
    FILELIST="$OUTPUT_DIR/genome_list.txt"
    > "$FILELIST"
    
    GENOME_COUNT=0
    for genome in "$GENOME_DIR"/*.{fna,fasta,fa}; do
        if [ -f "$genome" ]; then
            echo "$genome" >> "$FILELIST"
            GENOME_COUNT=$((GENOME_COUNT + 1))
            if [ $MAX_GENOMES -gt 0 ] && [ $GENOME_COUNT -ge $MAX_GENOMES ]; then
                break
            fi
        fi
    done
    
    echo "Файлов для обработки: $GENOME_COUNT"
    
    # Запуск AIndex
    echo "Построение индекса AIndex..."
    ./aindex build -k $K -t $THREADS -i "$FILELIST" -o "$OUTPUT_DIR/genome.aidx"
    
    # Анализ нулломеров
    echo "Анализ нулломеров..."
    ./aindex nullomers -i "$OUTPUT_DIR/genome.aidx" -o "$OUTPUT_DIR/aindex_nullomers.txt"
    
    echo -e "${GREEN}✓ AIndex анализ завершен${NC}"
}

# Функция для визуализации результатов
visualize_results() {
    echo -e "${YELLOW}Создание визуализаций...${NC}"
    
    # Python скрипт для визуализации будет создан отдельно
    if [ -f "visualize_nullomers.py" ]; then
        python3 visualize_nullomers.py "$OUTPUT_DIR" --k $K
    fi
}

# Основной процесс
echo -e "${GREEN}=== Пайплайн анализа нулломеров ===${NC}"
echo -e "${BLUE}Метод:${NC} $METHOD"
echo -e "${BLUE}k-мер:${NC} $K"
echo -e "${BLUE}Потоков:${NC} $THREADS"
echo -e "${BLUE}Выходная директория:${NC} $OUTPUT_DIR"
echo ""

# Компиляция программ
if [ "$SKIP_BUILD" = false ]; then
    compile_programs
fi

# Запуск выбранного метода
START_TIME=$(date +%s)

case $METHOD in
    cpp)
        run_cpp_method
        ;;
    jellyfish)
        run_jellyfish_method
        ;;
    aindex)
        run_aindex_method
        ;;
    *)
        echo -e "${RED}Неизвестный метод: $METHOD${NC}"
        exit 1
        ;;
esac

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo -e "${GREEN}=== Анализ завершен ===${NC}"
echo -e "Время выполнения: $(($ELAPSED / 60)) мин $(($ELAPSED % 60)) сек"
echo -e "Результаты сохранены в: ${BLUE}$OUTPUT_DIR${NC}"

# Создание сводного отчета
echo -e "\n${YELLOW}Создание сводного отчета...${NC}"
{
    echo "# Отчет по анализу нулломеров"
    echo "Дата: $(date)"
    echo "Метод: $METHOD"
    echo "k-мер: $K"
    echo "Геномов обработано: $GENOME_COUNT"
    echo "Время выполнения: $(($ELAPSED / 60)) мин $(($ELAPSED % 60)) сек"
    echo ""
    echo "## Файлы результатов:"
    ls -la "$OUTPUT_DIR"
} > "$OUTPUT_DIR/report.txt"

# Визуализация результатов
visualize_results

echo -e "${GREEN}✓ Все задачи выполнены${NC}"
