
# 1)Write a function that returns a float corresponding to the volume of a sphere: get_sphere_volume(radius) (formula=(4/3)πr3)

import math
r = float(input("Sphere radius: "))
sphere_volume = (4/3)*math.pi*r**3
print("The sphere volume is:", sphere_volume)


# 2) Write a function that calculates and returns an integer corresponding to the factorial of an integer (n):

    # a) Using recursivity: recursive_factorial(n)

n = int(input("Write an integer to calculate its factorial: ")) #The user is prompted to enter an integer n

def recursive_factorial(n):
    if n < 0:
        return "Negative numbers are not accepted" #The factorial function can be applied to all real numbers, except negative integers
    elif n == 0: #the factorial of 0 (an empty product) is 1
        return 1
    else: 
        return n*recursive_factorial(n-1)

result = recursive_factorial(n) #result variable accumulates the product of the numbers from 1 to n
print(f"The factorial of {n} is:", result)

    # b) Without using recursivity: factorial(n)

    # A for loop is used to iterate from 1 to n. 
    # For each iteration, the result variable previously initialized will be the same variable multiplied by the current iteration.

def factorial(n):
    result = 1
    for i in range(1, n + 1):
        result *= i
    return result

n = int(input("Write an integer to calculate its factorial: "))
print(f"The factorial of {n} is:", factorial(n))

# 3) Write a function for counting up numbers from 0 to n, showing the count up in the screen. If parameter odd is set to True, prints only odd numbers

    # a) Using recursivity: recursive_count_up(n, odd)
import math
def recursive_count_up(i, n, odd):
    if i > n:
        return
    if odd:
        if i % 2 != 0: #If the variable odd takes the value true, there is an "if" verifying if the remainder of the current value "i" divided by 2 is not equal to 0, therefore the result is displayed
            print(i)
    else:
        print(i)
    recursive_count_up(i + 1, n, odd) #the value of the iteration "i" sums 1

n = int(input("Write a number: "))
i = 0
odd = True #If the odd variable is assigned the value of true, only the odd values of the counter from 0 to n will appear on the screen
recursive_count_up(i, n, odd)

    # b) Without using recursivity: count_up(n,odd)

import math
def count_up(n, odd=False):
    for i in range(n + 1):
        if odd and i % 2 != 0:
                print(i)
        elif not odd:
            print(i)

n = int(input("Write a number: "))
odd_choise = input("Counting only odd numbers? (true/false): ").lower() #lowercase letters
odd = odd_choise == 'true' #boolean value
count_up(n, odd)

# 4) Find and solve the bugs in the following function: 
# def get_final_price(discount_percentage=10, price): """Return the final price after applying the discount percentage """
#return ((price + price) * percentage) / 100

def get_final_price(price, discount_percentage=10):
    
    """Return the final price after applying the discount percentage """
    return price - ((price * discount_percentage) / 100)

n = int(input("Write a price: "))
print(f"The final price is:", get_final_price(n))