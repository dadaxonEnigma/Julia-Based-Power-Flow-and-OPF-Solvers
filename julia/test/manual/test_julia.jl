# Первый тест Julia
println("="^60)
println("HELLO FROM JULIA!")
println("="^60)

# Базовый синтаксис Julia
println("\n1. BASIC OPERATIONS:")
a = 5
b = 10
println("a = $a, b = $b")
println("a + b = $(a + b)")

# Массивы
println("\n2. ARRAYS:")
arr = [1, 2, 3, 4, 5]
println("Array: $arr")
println("Sum: $(sum(arr))")

# Линейная алгебра
using LinearAlgebra

println("\n3. LINEAR ALGEBRA:")
A = [2.0 -1.0; -1.0 2.0]
b_vec = [1.0, 1.0]
x = A \ b_vec  # Решение системы Ax = b
println("Matrix A:")
println(A)
println("\nVector b: $b_vec")
println("Solution x = A\\b: $x")
println("Verification A*x: $(A*x)")

println("\n" * "="^60)
println("JULIA TEST COMPLETED!")
println("="^60)