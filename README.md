# MO11 PROJEKT

[mo_polecenie.pdf](https://github.com/Kuboss949/MO11/files/11824493/mo_polecenie.pdf)



UWAGA:
Jak wygląda oddawanie projektu?
Ogólnie to Pan dr bierze po jednej osobie, przegląda sprawozdanie i jeśli wszystko jest w miarę git, może coś lekko zapytać, ale jakoś szczególnie dobrze odpowiadać nie trzeba. 
Jeśli jest mało rzeczy do których się może doczepić to daje zwykle 4/4.5. Ogólnie oddawanie jest dość przyjemne. Nie przegląda jakoś szczególnie kodu, chyba, że mu się coś bardzo niepodoba.

0.5 obcina za brak wyjaśnienia czemu błąd zmienia się wraz z t.

WYJAŚNIENIE czemu:
W moim wypadku funkcja w warunku początkowym jest nieciągła, co powoduje, że nie ma pochodnej, dlatego błędy na początku są tak duże, bo przybliżamy pochodną, która jest nieokreślona. 
Ogólnym problemem tutaj jest błąd obcięcia dla naszego równania. Zależy on od kolejnych pochocnych naszej funkcji (jak przy przybliżeniach różnicowych), które mogą zachowywać się różnie, np:

Może być tak, że funkcja na początku jest "wybrzuszona", a potem maleje

![obraz](https://github.com/Kuboss949/MO11/assets/101654879/b7249935-330b-48d3-8e8b-aa36d048beaa)

Mniej więcej coś takiego. Wtedy funkcja "dąży do prostej", a dalsze pochodne do 0, więc błąd też dąży do zera

Może być też tak, że funkcja dąży do pewnego stanu ustalonego

![obraz](https://github.com/Kuboss949/MO11/assets/101654879/d0e17ab0-ca02-4228-8810-19eb2b0341ca)

czerwona linia - stan ustalony, żółta - początkowy

Wtedy pochodne będą stałe, przez co błąd też w pewnym momencie stanie się stały.

CZEGO UNIKAĆ

Słowa stabilizacja (nie lubi tego w sprawozdaniu) :)

Nie używania biblioteki, którą udostępnia

