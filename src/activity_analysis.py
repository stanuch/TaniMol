import numpy as np


def within_cluster_activity_distributions(clusters: dict, pic50: np.ndarray):
    results = {}
    for centroid, members in clusters.items():
        values = pic50[list(members)]
        results[centroid] = {
            "n": len(values),
            "mean": np.mean(values),
            "median": np.median(values),
            "std": np.std(values),
            "min": np.min(values),
            "max": np.max(values),
        }
    return results


def activity_cliffs():
    raise NotImplementedError


# NOTE:
# Dobra, tutaj trzeba chyba trochę bardziej pogłówkować tak mi się wydaje
# Trzeba znaleźć pary cząsteczek które mają wysokie similarity Tanimoto (np. 0.8)
# ale mające duże różnice w pIC50 (np. > 2).
# Czyli - trzeba sprawdzić w similarity matrix wartości > 0.8, następnie zobaczyć jakie cząsteczki tworzą tę parę.
# Później trzeba wczytać cleaned_activities.csv i dla tych dwóch cząsteczek sprawdzić pIC50 values.
# Jeśli różnica jest > 2, to znaczy, że jest activity cliff i trochę nie wiem co to znaczy. W sensie czy jakaś konkretna
# zmiana w cząstecze to spowodowała, czy da się to jakoś zwizalizować. Nawet jeśli dostanę dane, to trochę nie wiem
# co można z nimi dalej zrobić. Wydaje mi się to dobrym podejściem, aczkolwiek sprawdzanie od początku
# całego matrixa może być czasochłonne i zasobochłonne. Więc w sumie nie wiem czy to podejście jest okej.
# Pierwsze jakie mi przyszło do głowy w sumie.


def sali():
    raise NotImplementedError


# NOTE:
# Tutaj tylko zamiast Tanimoto, mamy dystans, ale to w sumie to samo, tylko że odwrócone.
# To jest też chyba łatwiejsze do sprawdzenia i zwizualiowania tak mi się wydaje. Bo nie muszę pokazywać
# Które to są dokładnie cząsteczki i z nich robić wykresu jakiegoś, tylko po prostu biorę sobie
# każdą cząsteczke i obliczam dla pary (z inną cząsteczką) z nią współczynnik SALI. Wtedy można go wrzucić na histogram jakiś
# i wtedy widać, czy dużo jest takich cząsteczek w których jedna zmiana powoduje dużą różnicę
# Bo chyba generalnie o to w tym chodzi. Że sprawdzamy jak jakaś drobna (dla Tanimoto) zmiana, powoduje
# duże różnice w aktwności.


def similarity_activity_correlation():
    raise NotImplementedError


# NOTE:
# Mam dwie macierze NxN: similarity_matrix (Tanimoto) i delta_pic50 (|ΔpIC50|).
# Obie są symetryczne, więc biorę tylko górny trójkąt (żeby nie liczyć par dwa razy).
# To daje mi dwa równoległe wektory:
#   - wektor_1: wartości Tanimoto dla każdej pary (i,j)
#   - wektor_2: |ΔpIC50| dla tej samej pary (i,j)
# Wrzucam oba wektory do Spearmana i dostaję jeden współczynnik korelacji.
# Jeśli korelacja jest ujemna → im wyższe Tanimoto (bardziej podobne struktury),
# tym niższe |ΔpIC50| (bardziej zbliżona aktywność). To potwierdza hipotezę SAR.
# Jeśli korelacja jest bliska zeru → struktura nie przewiduje aktywności,
# co oznacza że SAR nie działa dla tego zbioru danych.


# NOTE:
# Generalnie zastanawiam się jeszcze czy to nie jest trochę redundant i marnowanie to co robię
# Mam wrażenie, jakby dało się to zrobić wydajniej a nie co chwilę przelatywać cały matrix od początku
# i robić kurwa "hope for the best". Nie jestem pewien czy te rozwiązania są okej.
