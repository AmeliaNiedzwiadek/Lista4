class Wielomian:
    def __init__(self, wspolczynniki):
        """
        Konstruktor klasy Wielomian.

        @param wspolczynniki: lista współczynników, gdzie i-ty element odpowiada współczynnikowi przy x^i
        @throws ValueError: jeśli lista współczynników jest pusta
        """
        if not wspolczynniki:
            raise ValueError("Lista współczynników nie może być pusta.")
        while len(wspolczynniki) > 1 and wspolczynniki[-1] == 0:
            wspolczynniki.pop()
        self.wsp = wspolczynniki

    def stopien(self):
        """
        Zwraca stopień wielomianu.

        @return: int - stopień wielomianu
        """
        return len(self.wsp) - 1

    def __str__(self):
        """
        Zwraca reprezentację tekstową wielomianu.

        @return: str - postać tekstowa wielomianu
        """
        czesci = []
        for i in range(len(self.wsp) - 1, -1, -1):
            wsp = self.wsp[i]
            if wsp == 0:
                continue
            znak = "+" if wsp > 0 else "-"
            wsp_abs = abs(wsp)
            if i == 0:
                czesci.append(f"{znak} {wsp_abs}")
            elif i == 1:
                czesci.append(f"{znak} {wsp_abs}x")
            else:
                czesci.append(f"{znak} {wsp_abs}x^{i}")
        if not czesci:
            return "W(x) = 0"
        czesci[0] = czesci[0].replace("+", "")  # pierwszy znak bez +
        return "W(x) = " + " ".join(czesci)

    def __call__(self, x):
        """
        Oblicza wartość wielomianu dla danej wartości x.

        @param x: float - argument
        @return: float - wartość wielomianu w punkcie x
        """
        return sum(self.wsp[i] * (x ** i) for i in range(len(self.wsp)))

    def __add__(self, other):
        """
        Dodaje dwa wielomiany.

        @param other: Wielomian
        @return: Wielomian - suma dwóch wielomianów
        """
        maks = max(len(self.wsp), len(other.wsp))
        wynik = []
        for i in range(maks):
            a = self.wsp[i] if i < len(self.wsp) else 0
            b = other.wsp[i] if i < len(other.wsp) else 0
            wynik.append(a + b)
        return Wielomian(wynik)

    def __sub__(self, other):
        """
        Odejmuje wielomian other od bieżącego.

        @param other: Wielomian
        @return: Wielomian - różnica wielomianów
        """
        maks = max(len(self.wsp), len(other.wsp))
        wynik = []
        for i in range(maks):
            a = self.wsp[i] if i < len(self.wsp) else 0
            b = other.wsp[i] if i < len(other.wsp) else 0
            wynik.append(a - b)
        return Wielomian(wynik)

    def __mul__(self, other):
        """
        Mnoży dwa wielomiany.

        @param other: Wielomian
        @return: Wielomian - iloczyn wielomianów
        """
        wynik = [0] * (len(self.wsp) + len(other.wsp) - 1)
        for i in range(len(self.wsp)):
            for j in range(len(other.wsp)):
                wynik[i + j] += self.wsp[i] * other.wsp[j]
        return Wielomian(wynik)



# TESTY 
def test(opis, kod):
    print(f"\n{opis}")
    try:
        kod()
    except Exception as e:
        print(" Błąd:", e)

print(" Testy KLASY 'Wielomian'")

test("1. Tworzenie wielomianu [1, 2, 3]:", lambda: print(Wielomian([1, 2, 3])))

test("2. Tworzenie pustego wielomianu:",
     lambda: print(Wielomian([])))

test("3. Skracanie zer z końca: [5, 0, 0] ⇒ [5]", lambda: print(Wielomian([5, 0, 0])))

test("4. Stopień wielomianu [3, 0, 2]:",
     lambda: print(f"{Wielomian([3, 0, 2])} | stopień = {Wielomian([3, 0, 2]).stopien()}"))

test("5. Obliczanie wartości w punkcie x = 1 dla [1, 2, 1]:",
     lambda: print(f"{Wielomian([1, 2, 1])}(1) = {Wielomian([1,2,1])(1)}"))

test("6. Dodawanie [1, 2] + [3, 4, 5]:",
     lambda: print(Wielomian([1, 2]), "+", Wielomian([3, 4, 5]), "=", Wielomian([1, 2]) + Wielomian([3, 4, 5])))

test("7. Odejmowanie [4, 5] - [1, 2]:",
     lambda: print(Wielomian([4, 5]), "-", Wielomian([1, 2]), "=", Wielomian([4, 5]) - Wielomian([1, 2])))

test("8. Mnożenie [1, 1] * [1, 1]:",
     lambda: print(Wielomian([1, 1]), "*", Wielomian([1, 1]), "=", Wielomian([1, 1]) * Wielomian([1, 1])))

test("9. Wielomian zerowy [0, 0, 0]:", lambda: print(Wielomian([0, 0, 0])))

test("10. Wyświetlanie znaku minus: [1, 0, -3]:", lambda: print(Wielomian([1, 0, -3])))

test("11. Dodanie do zera: [0] + [2,1]:", lambda: print(Wielomian([0]) + Wielomian([2, 1])))

test("12. Różnica dwóch takich samych: [1,2,3] - [1,2,3]:",
     lambda: print(Wielomian([1, 2, 3]) - Wielomian([1, 2, 3])))

test("13. Mnożenie z zerem: [1,2] * [0]:", lambda: print(Wielomian([1, 2]) * Wielomian([0])))

test("14. Wielomian stopnia 0: [5]", lambda: print(Wielomian([5])))

test("15. Wielomian ujemny: [-2, -1, -3]:", lambda: print(Wielomian([-2, -1, -3])))

test("16. Złożona operacja: ([1,2] + [0,1]) * [1,0,1]:",
     lambda: print((Wielomian([1, 2]) + Wielomian([0, 1])) * Wielomian([1, 0, 1])))

test("17. Długi wielomian: [1,1,1,1,1,1]:", lambda: print(Wielomian([1, 1, 1, 1, 1, 1])))

test("18. Oblicz wartość wielomianu x^3 + 2x^2 - x + 1 w x=2:",
     lambda: print(Wielomian([1, -1, 2, 1])(2)))

test("19. Mnożenie dwóch różnych długości: [1,2] * [1]:",
     lambda: print(Wielomian([1, 2]) * Wielomian([1])))

#Źródła z których korzystałam tworząc ten kod to:
#https://docs.python.org/3/reference/datamodel.html#basic-customization
#https://www.oracle.com/technical-resources/articles/java/javadoc-tool.html
#https://docs.python.org/3/reference/expressions.html#lambda
#ChatGPT
#https://realpython.com/python-lambda/
#https://www.learnpython.org/en/Exception_Handling
#https://docs.python.org/3/
#https://docs.python.org/3/library/dataclasses.html
#https://github.com/AmeliaNiedzwiadek/lista-2/blob/main/zadanie1.kt


