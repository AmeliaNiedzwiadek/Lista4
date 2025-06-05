class BioSequence:
    """
    Klasa bazowa reprezentująca ogólną sekwencję biologiczną (np. DNA, RNA, białka).

    Atrybuty:
        identifier (str): Identyfikator sekwencji (np. nazwa lub opis).
        data (str): Sekwencja znaków (np. nukleotydów lub aminokwasów), przechowywana wielkimi literami.
        valid_chars (set): Zbiór dozwolonych znaków w sekwencji.

    Wyjątki:
        ValueError: Gdy sekwencja zawiera niedozwolone znaki.
    """

    def __init__(self, identifier: str, data: str, valid_chars: set):
        """
        Inicjalizuje obiekt sekwencji biologicznej i sprawdza poprawność znaków.

        Args:
            identifier (str): Identyfikator sekwencji.
            data (str): Surowa sekwencja.
            valid_chars (set): Dozwolone znaki.

        Raises:
            ValueError: Jeśli `data` zawiera znaki spoza `valid_chars`.
        """
        self.identifier = identifier
        self.data = data.upper()
        self.valid_chars = valid_chars
        if not all(char in valid_chars for char in self.data):
            raise ValueError(f"Nieprawidłowe znaki w sekwencji: {self.data}")

    def length(self):
        """
        Zwraca długość sekwencji.

        Returns:
            int: Liczba znaków w sekwencji.
        """
        return len(self.data)

    def __str__(self):
        """
        Zwraca reprezentację sekwencji w formacie FASTA.

        Returns:
            str: Reprezentacja tekstowa.
        """
        return f">{self.identifier}\n{self.data}"

    def mutate(self, position: int, value: str):
        """
        Zmienia znak w sekwencji na podanej pozycji.

        Args:
            position (int): Indeks (0-indeksowany), który ma zostać zmodyfikowany.
            value (str): Nowy znak (musi być dozwolony).

        Raises:
            IndexError: Jeśli `position` jest poza zakresem sekwencji.
            ValueError: Jeśli `value` nie należy do `valid_chars`.
        """
        if position < 0 or position >= len(self.data):
            raise IndexError("Pozycja poza zakresem.")
        if value not in self.valid_chars:
            raise ValueError(f"Znak {value} nie jest dozwolony.")
        self.data = self.data[:position] + value + self.data[position+1:]

    def find_motif(self, motif: str) -> int:
        """
        Wyszukuje motyw (podciąg) w sekwencji.

        Args:
            motif (str): Podciąg do wyszukania.

        Returns:
            int: Indeks pierwszego wystąpienia motywu lub -1, jeśli nie znaleziono.
        """
        return self.data.find(motif)


class DNASequence(BioSequence):
    """
    Klasa reprezentująca sekwencję DNA. Dziedziczy z BioSequence i definiuje znaki 'A', 'T', 'G', 'C'.
    """

    def __init__(self, identifier: str, data: str):
        """
        Inicjalizuje obiekt DNA z odpowiednim zestawem dozwolonych znaków.

        Args:
            identifier (str): Identyfikator sekwencji.
            data (str): Sekwencja DNA.
        """
        super().__init__(identifier, data, {'A', 'T', 'G', 'C'})

    def complement(self) -> str:
        """
        Zwraca komplementarną (i odwróconą) sekwencję DNA.

        Returns:
            str: Komplementarna sekwencja DNA.
        """
        complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        comp = ''.join(complement_map[base] for base in reversed(self.data))
        return comp

    def transcribe(self):
        """
        Inicjalizuje obiekt RNA na podstawie komplementarnej, odwróconej sekwencji DNA.

        Returns:
            RNASequence: Nowa sekwencja RNA.
        """
        complement_map = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}
        rna_data = ''.join(complement_map[base] for base in reversed(self.data))
        return RNASequence(self.identifier, rna_data)


class RNASequence(BioSequence):
    """
    Klasa reprezentująca sekwencję RNA. Używa znaków 'A', 'U', 'G', 'C'.
    """

    def __init__(self, identifier: str, data: str):
        """
        Inicjalizuje sekwencję RNA z odpowiednim zestawem znaków.

        Args:
            identifier (str): Identyfikator sekwencji.
            data (str): Sekwencja RNA.
        """
        super().__init__(identifier, data, {'A', 'U', 'G', 'C'})

    def complement(self) -> str:
        """
        Zwraca komplementarną (i odwróconą) sekwencję RNA.

        Returns:
            str: Komplementarna sekwencja RNA.
        """
        complement_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
        comp = ''.join(complement_map[base] for base in reversed(self.data))
        return comp

    def transcribe(self):
        """
        Tłumaczy sekwencję RNA na sekwencję aminokwasową (białko), do pierwszego kodonu stop.

        Returns:
            ProteinSequence: Obiekt reprezentujący białko.
        """
        codon_table = {
            'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
            'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
            'UAU': 'Y', 'UAC': 'Y', 'UGU': 'C', 'UGC': 'C',
            'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
            'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
            'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
            'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
            'UAA': '*', 'UAG': '*', 'UGA': '*'
        }

        protein = ''
        for i in range(0, len(self.data) - 2, 3):
            codon = self.data[i:i+3]
            amino_acid = codon_table.get(codon, 'X')
            if amino_acid == '*':
                break
            protein += amino_acid
        return ProteinSequence(self.identifier, protein)


class ProteinSequence(BioSequence):
    """
    Klasa reprezentująca sekwencję białkową. Używa jednoznakowych oznaczeń aminokwasów.
    """

    def __init__(self, identifier: str, data: str):
        """
        Inicjalizuje sekwencję białkową.

        Args:
            identifier (str): Identyfikator sekwencji.
            data (str): Sekwencja aminokwasów (np. 'MKT...').
        """
        allowed = set("ACDEFGHIKLMNPQRSTVWY")
        super().__init__(identifier, data, allowed)


# ===================== TESTY ======================

def assert_equals(expected, actual, test_name):
    result = "OK" if expected == actual else f"FAIL — Expected: {expected}, Actual: {actual}"
    print(f"{test_name}: {result}")

def main():
    dna = DNASequence("DNA1", "ATGC")
    test_basic(dna)
    test_mutations(dna)
    test_motifs(dna)
    test_transcription_and_complement(dna)

def test_basic(dna):
    print("\n=== Test podstawowy ===")
    print(dna)
    print("Długość:", dna.length)

def test_mutations(dna):
    print("\n=== Testowanie mutacji ===")

    print("-> Mutacja pozycji 2 na 'T'")
    try:
        dna.mutate(2, 'T')
        print("Po mutacji:", dna)
        assert_equals("ATTC", dna.data, "Mutacja poprawna")
    except Exception as e:
        print(" Nieoczekiwany błąd przy mutacji pozycji 2:", e)

    print("-> Próba mutacji na pozycji 10 ")
    try:
        dna.mutate(10, 'A')
        print(" Błąd! Oczekiwano wyjątku, ale mutacja się powiodła.")
    
    except Exception as e:
        print(f" Niespodziewany typ wyjątku: {e}")

    print("-> Próba mutacji pozycji 0 na niedozwolony znak 'X'")
    try:
        dna.mutate(0, 'X')
        print(" Błąd! Oczekiwano wyjątku, ale mutacja się powiodła.")
    except ValueError as e:
        print(f" {e}")
    except Exception as e:
        print(f" Niespodziewany typ wyjątku: {e}")


def test_motifs(dna):
    print("\n=== Testowanie motywów ===")
    pos = dna.find_motif("TT")
    print("Motyw 'TT' na pozycji:", pos)
    assert_equals(1, pos, "Motyw TT")

    no_pos = dna.find_motif("GGG")
    print("Motyw 'GGG' na pozycji:", no_pos)
    assert_equals(-1, no_pos, "Brak motywu")

def test_transcription_and_complement(dna):
    print("\n=== Transkrypcja i komplement ===")
    comp = dna.complement()
    print("Komplement DNA:", comp)
    assert_equals("GAAT", comp, "Komplement DNA (odwrócony)")

    rna = dna.transcribe()
    print(f"RNA (z {dna.data}): {rna.data}")
    assert_equals("GAAU", rna.data, "Transkrypcja RNA")

    rna_comp = rna.complement()
    print("Komplement RNA:", rna_comp)
    assert_equals("AUUC", rna_comp, "Komplement RNA")

    protein = rna.transcribe()
    print("Białko:", protein)
    print("Sekwencja białka:", protein.data)

if __name__ == "__main__":
    main()

#Podczas tworzenia kodu korzystałam z następujących źródeł:
#https://docs.python.org/3/tutorial/classes.html
#https://docs.python.org/3/library/functions.html#super
#https://docs.python.org/3/library/functions.html#all
#https://docs.python.org/3/library/stdtypes.html#str.find
#https://docs.python.org/3/library/stdtypes.html#text-sequence-type-str
#https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#https://realpython.com/python3-object-oriented-programming/
#ChatGPT
#https://docs.python.org/3/tutorial/inputoutput.html#formatted-string-literals
#https://docs.python.org/3/tutorial/errors.html
#https://docs.python.org/3/reference/compound_stmts.html#the-try-statement
#https://realpython.com/python-testing/#writing-basic-tests-with-assert
#https://realpython.com/python-main-function/
#https://docs.python.org/3/reference/expressions.html#comparisons
#https://github.com/AmeliaNiedzwiadek/lista-2
