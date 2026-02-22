#!/bin/bash

# === KONFIGURATION ===
SOURCE_DIR="$1"
TARGET_DIR="$2"

if [ -z "$SOURCE_DIR" ] || [ -z "$TARGET_DIR" ]; then
    echo "Usage: $0 <source_dir> <target_dir>"
    exit 1
fi

REQUIRED_FILES=("E_field.vtu" "metadata.json")
VISU_DIR="${TARGET_DIR}visu"
mkdir -p "$VISU_DIR" &>/dev/null


echo "========================================================================"
echo "MaxwellScattererAndAntennaNGSolve Clean-Up Dataset"
echo "========================================================================"
echo "Source Path: ${SOURCE_DIR}"
echo "Target Path: ${TARGET_DIR}"
echo ""
echo "Required Files:"
for file in "${REQUIRED_FILES[@]}"; do
    echo "  - ${file}"
done
echo ""

# =====================

# Unterordner zählen
mapfile -t DIRS < <(find "$SOURCE_DIR" -mindepth 1 -maxdepth 1 -type d)
TOTAL=${#DIRS[@]}

if [ "$TOTAL" -eq 0 ]; then
    echo "Keine Unterordner gefunden."
    exit 0
fi

# === Farben ===
RESET=$'\033[0m'
RED=$'\033[31m'
YELLOW=$'\033[33m'
GREEN=$'\033[32m'
CYAN=$'\033[36m'
WHITE=$'\033[37m'

# =====================

current=0
fails=0
htmls=0
BLOCK_SIZE=3
NUM_BLOCKS=$(( (TOTAL + BLOCK_SIZE - 1) / BLOCK_SIZE ))

update_bar_width() {
    TERM_WIDTH=$(tput cols)
    RESERVED=25
    BAR_WIDTH=$((TERM_WIDTH - RESERVED))
    [ "$BAR_WIDTH" -lt 10 ] && BAR_WIDTH=10
}
update_bar_width
trap update_bar_width SIGWINCH

tput civis

# === Platz für Progressbar reservieren ===
PROGRESS_LINES=$((NUM_BLOCKS + 6))
printf "\n%.0s" $(seq 1 $PROGRESS_LINES)

draw_progress() {
    folder="$1"
    percent_total=$((current * 100 / TOTAL))

    # Farben für Gesamt
    if [ "$percent_total" -lt 50 ]; then total_color=$RED
    elif [ "$percent_total" -lt 90 ]; then total_color=$YELLOW
    else total_color=$GREEN
    fi

    printf "\033[%dA" "$PROGRESS_LINES"

    # Rahmen oben
    printf "\033[2K\r${CYAN}╔"
    printf '═%.0s' $(seq 1 $((TERM_WIDTH-2)))
    printf "╗${RESET}\n"

    # Hilfsfunktion: saubere Zeilen
    print_row() {
        visible="$1"
        colored="$2"
        padding=$((TERM_WIDTH - ${#visible} - 3))
        [ "$padding" -lt 0 ] && padding=0
        printf "\033[2K\r${CYAN}║${RESET} %b%*s${CYAN}║${RESET}\n" "$colored" "$padding" ""
    }

    # Gesamtzeile mit HTML/Fails
    visible_total="Gesamt: $current / $TOTAL (${percent_total}%) | HTML: $htmls | Fails: $fails"
    total_line="Gesamt: $current / $TOTAL (${total_color}${percent_total}%${RESET}) | HTML:${GREEN} $htmls${RESET} | Fails:${RED} $fails${RESET}"
    print_row "$visible_total" "$total_line"

    # Leerzeile
    print_row "" ""

    # Blocks
    for ((b=0; b<NUM_BLOCKS; b++)); do
        start=$((b * BLOCK_SIZE + 1))
        end=$((start + BLOCK_SIZE - 1))
        [ "$end" -gt "$TOTAL" ] && end="$TOTAL"
        block_total=$((end - start + 1))

        if [ "$current" -lt "$start" ]; then block_pos=0
        elif [ "$current" -gt "$end" ]; then block_pos=$block_total
        else block_pos=$((current - start + 1))
        fi

        percent_block=$((block_pos * 100 / block_total))
        if [ "$percent_block" -lt 50 ]; then color=$RED
        elif [ "$percent_block" -lt 90 ]; then color=$YELLOW
        else color=$GREEN
        fi

        filled=$((block_pos * BAR_WIDTH / block_total))
        empty=$((BAR_WIDTH - filled))
        bar_filled=$(printf "%${filled}s" | tr ' ' '#')
        bar_empty=$(printf "%${empty}s")

        visible_block="Block $((b+1)) ($start–$end) [${bar_filled}${bar_empty}] ${percent_block}%"
        colored_block="Block $((b+1)) ($start–$end) [${color}${bar_filled}${RESET}${bar_empty}] ${color}${percent_block}%${RESET}"
        print_row "$visible_block" "$colored_block"
    done

    # Leerzeile
    print_row "" ""

    # Aktueller Ordner
    current_line="Aktueller Ordner: $folder"
    print_row "$current_line" "$current_line"

    # Rahmen unten
    printf "\033[2K\r${CYAN}╚"
    printf '═%.0s' $(seq 1 $((TERM_WIDTH-2)))
    printf "╝${RESET}\n"
}

# === Hauptschleife ===
for dir in "${DIRS[@]}"; do
    [ -d "$dir" ] || continue
    ((current++))

    folder_name=$(basename "$dir")
    job_id="${folder_name#job_}"

    # Fortschrittsbalken zuerst
    draw_progress "$folder_name"

    # HTML-Datei verschieben
    html_file="$dir/field_visualization.html"
    if [ -f "$html_file" ]; then
        new_file="field_visualization_${job_id}.html"
        ((htmls++))
        mv "$html_file" "$VISU_DIR/$new_file" &>/dev/null
    fi

    # Überprüfe erforderliche Dateien
    all_files_exist=true
    for file in "${REQUIRED_FILES[@]}"; do
        if [ ! -f "$dir/$file" ]; then
            all_files_exist=false
            ((fails++))
            break
        fi
    done

    # Dateien verschieben, wenn alle vorhanden
    if [ "$all_files_exist" = true ]; then
        target_folder="$TARGET_DIR$folder_name"
        mkdir -p "$target_folder" &>/dev/null
        for file in "${REQUIRED_FILES[@]}"; do
            mv "$dir/$file" "$target_folder/$file" &>/dev/null
        done
    fi
done

# Cursor wieder unterhalb der Progressbar
printf "\033[%dB" "$PROGRESS_LINES"
tput cnorm

echo -e "\nFertig!"

# #!/bin/bash

# # === KONFIGURATION ===
# SOURCE_DIR="$1"
# TARGET_DIR="$2"

# if [ -z "$SOURCE_DIR" ] || [ -z "$TARGET_DIR" ]; then
#     echo "Usage: $0 <source_dir> <target_dir>"
#     exit 1
# fi

# # Dateien, die vorhanden sein müssen
# REQUIRED_FILES=("E_field.vtu" "metadata.json")
# # REQUIRED_FILES=("E_field.vtk" "E_field.vtu" "metadata.json")

# echo "========================================================================"
# echo "MaxwellScattererAndAntennaNGSolve Clean-Up Dataset"
# echo "========================================================================"
# echo "Source Path: ${SOURCE_DIR}"
# echo "Target Path: ${TARGET_DIR}"
# echo ""
# echo "Required Files:"
# for file in "${REQUIRED_FILES[@]}"; do
#     echo "  - ${file}"
# done
# echo ""

# VISU_DIR="${TARGET_DIR}visu"
# mkdir -p "$VISU_DIR"

# # =====================

# # Unterordner zählen
# mapfile -t DIRS < <(find "$SOURCE_DIR" -mindepth 1 -maxdepth 1 -type d)
# TOTAL=${#DIRS[@]}

# if [ "$TOTAL" -eq 0 ]; then
#     echo "Keine Unterordner gefunden."
#     exit 0
# fi

# # === Farben ===
# RESET=$'\033[0m'
# BOLD=$'\033[1m'

# RED=$'\033[31m'
# YELLOW=$'\033[33m'
# GREEN=$'\033[32m'
# CYAN=$'\033[36m'
# WHITE=$'\033[37m'

# # =====================

# current=0
# fails=0
# htmls=0

# BLOCK_SIZE=3
# NUM_BLOCKS=$(( (TOTAL + BLOCK_SIZE - 1) / BLOCK_SIZE ))

# update_bar_width() {
#     TERM_WIDTH=$(tput cols)
#     RESERVED=25
#     BAR_WIDTH=$((TERM_WIDTH - RESERVED))
#     [ "$BAR_WIDTH" -lt 10 ] && BAR_WIDTH=10
# }
# update_bar_width
# trap update_bar_width SIGWINCH

# tput civis

# # Reserven für Progressbar-Frame
# PROGRESS_LINES=$((NUM_BLOCKS + 6))
# # Cursor nach unten bewegen, um Platz für Progressbar zu reservieren
# printf "\n%.0s" $(seq 1 $PROGRESS_LINES)

# LINES=$((NUM_BLOCKS + 6))

# for ((i=0;i<LINES;i++)); do
#     echo
# done

# draw_progress() {

#     percent_total=$((current * 100 / TOTAL))

#     if [ "$percent_total" -lt 50 ]; then
#         total_color=$RED
#     elif [ "$percent_total" -lt 90 ]; then
#         total_color=$YELLOW
#     else
#         total_color=$GREEN
#     fi

#     printf "\033[%dA" "$LINES"

#     # Rahmen oben
#     printf "\033[2K\r${CYAN}╔"
#     printf '═%.0s' $(seq 1 $((TERM_WIDTH-2)))
#     printf "╗${RESET}\n"

#     # Hilfsfunktion für saubere Zeilen
#     print_row() {
#         visible="$1"
#         colored="$2"

#         visible_length=${#visible}
#         padding=$((TERM_WIDTH - visible_length - 3))
#         [ "$padding" -lt 0 ] && padding=0

#         printf "\033[2K\r${CYAN}║${RESET} %b%*s${CYAN}║${RESET}\n" \
#             "$colored" "$padding" ""
#     }

#     # ===== Gesamt =====
#     visible_total="Gesamt: $current / $TOTAL (${percent_total}%) | HTML: $htmls | Fails: $fails"

#     total_line="Gesamt: $current / $TOTAL ("
#     total_line+="${total_color}${percent_total}%${RESET})"

#     # HTML grün, Fails rot
#     total_line+=" | HTML:${GREEN} $htmls${RESET}"
#     total_line+=" | Fails:${RED} $fails${RESET}"

#     print_row "$visible_total" "$total_line"

#     # Leerzeile
#     print_row "" ""

#     # ===== Blocks =====
#     for ((b=0; b<NUM_BLOCKS; b++)); do

#         block_start=$(( b * BLOCK_SIZE + 1 ))
#         block_end=$(( block_start + BLOCK_SIZE - 1 ))
#         [ "$block_end" -gt "$TOTAL" ] && block_end="$TOTAL"

#         block_total=$(( block_end - block_start + 1 ))

#         if [ "$current" -lt "$block_start" ]; then
#             block_pos=0
#         elif [ "$current" -gt "$block_end" ]; then
#             block_pos=$block_total
#         else
#             block_pos=$(( current - block_start + 1 ))
#         fi

#         percent_block=$(( block_pos * 100 / block_total ))

#         if [ "$percent_block" -lt 50 ]; then
#             color=$RED
#         elif [ "$percent_block" -lt 90 ]; then
#             color=$YELLOW
#         else
#             color=$GREEN
#         fi

#         filled=$(( block_pos * BAR_WIDTH / block_total ))
#         empty=$(( BAR_WIDTH - filled ))

#         bar_filled=$(printf "%${filled}s" | tr ' ' '#')
#         bar_empty=$(printf "%${empty}s")

#         visible_block="Block $((b+1)) ($block_start–$block_end) [${bar_filled}${bar_empty}] ${percent_block}%"

#         colored_block="Block $((b+1)) ($block_start–$block_end) [${color}${bar_filled}${RESET}${bar_empty}] ${color}${percent_block}%${RESET}"

#         print_row "$visible_block" "$colored_block"
#     done

#     # Leerzeile
#     print_row "" ""

#     # Aktueller Ordner
#     current_line="Aktueller Ordner: $1"
#     print_row "$current_line" "$current_line"

#     # Rahmen unten
#     printf "\033[2K\r${CYAN}╚"
#     printf '═%.0s' $(seq 1 $((TERM_WIDTH-2)))
#     printf "╝${RESET}\n"
# }

# for dir in "${DIRS[@]}"; do
#     # Nur wenn es wirklich ein Verzeichnis ist
#     [ -d "$dir" ] || continue

#     ((current++))

#     folder_name=$(basename "$dir")
#     job_id="${folder_name#job_}"

#     draw_progress "$(basename "$dir")"

#     # HTML-Datei verschieben
#     html_file="$dir/field_visualization.html"
#     if [ -f "$html_file" ]; then
#         new_file="field_visualization_${job_id}.html"
#         ((htmls++))
#         mv "$html_file" "$VISU_DIR/$new_file" 2>/dev/null  # Fehlerausgabe unterdrücken
#     fi

#     all_files_exist=true
#     for file in "${REQUIRED_FILES[@]}"; do
#         if [ ! -f "$dir/$file" ]; then
#             all_files_exist=false
#             ((fails++))
#             break
#         fi
#     done

#     # Nur wenn alle Dateien da sind
#     if [ "$all_files_exist" = true ]; then
#         target_folder="$TARGET_DIR$folder_name"
#         mkdir -p "$target_folder" 2>/dev/null
#         for file in "${REQUIRED_FILES[@]}"; do
#             mv "$dir/$file" "$target_folder/$file" 2>/dev/null
#         done
#     fi
# done

# # Cursor wieder korrekt unterhalb der Progressbar
# printf "\033[%dB" "$LINES"
# tput cnorm

# # Danach normale Ausgaben
# echo -e "\nFertig!"