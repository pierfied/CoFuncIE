T="$(date +%s)"

./par_convert_cart < full_data.txt

T="$(($(date +%s)-T))"
echo "Time in seconds: ${T}"