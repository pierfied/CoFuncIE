T="$(date +%s)"

./convert_cart < full_data.txt

T="$(($(date +%s)-T))"
echo "Time in seconds: ${T}"