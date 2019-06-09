import sys

def test(prvi, drugi):
	granica = 25
	prvi = int(prvi)
	drugi = int(drugi)
	rez = abs(prvi - drugi)
	if rez > granica: return False
	return True

def test3(prvi, drugi):
	granica = 3
	prvi = int(prvi)
	drugi = int(drugi)
	rez = abs(prvi - drugi)
	if rez > granica: return False
	return True


minimap_dat = open(sys.argv[1], 'r')
mapper_dat = open(sys.argv[2], 'r')
k_value = int(sys.argv[3])

content = minimap_dat.readlines()
content = [x.strip() for x in content]

velicinaMinimapa = len(content)

duplici = 0
mapa = {}

dupliciList = list()
dupliciSet = set()

brojNemaPocetakNiKraj = 0;

provjerenoMinimap = set()

for line in content:
	key = (line.split()[0], line.split()[5])

	provjerenoMinimap.add(key)

	if int(line.split()[2]) > k_value and (int(line.split()[1]) - int(line.split()[3])) > k_value:
		brojNemaPocetakNiKraj += 1

	if key not in mapa:
		mapa[key] = line
	else:
		pocetak1, kraj1 = line.split()[2], line.split()[3]
		otherLine = mapa[key]
		pocetak2, kraj2 = otherLine.split()[2], otherLine.split()[3]
		if (int(kraj1) - int(pocetak1)) > (int(kraj2) - int(pocetak2)):
			mapa[key] = line
		
		duplici += 1
		dupliciList.append(line)
		dupliciSet.add(key)

content = mapper_dat.readlines()
content = [x.strip() for x in content]

velicinaMappera = len(content)

brojUnutarGranice3 = 0
brojUnutarGranice25 = 0
brojIstih = 0
brojMinimapNema = 0
brojRazlicitStrand = 0
brojNemaNiceg = 0
error = 0

provjereno = set()

for line in content:
	ime, velicina, pocetakQ, krajQ, strand, REFime, REFvelicina, pocetakR, krajR, negleda1, negleda2, negleda3 = line.split()
	key = (ime, REFime)

	provjereno.add(key)
	
	barNesto = False

	if int(krajQ) > int(velicina):
		error += 1
		continue

	if key not in mapa:
		brojMinimapNema += 1
		continue

	ime2, velicina2, pocetakQ2, krajQ2, strand2, REFime2, REFvelicina2, pocetakR2, krajR2, negleda12, negleda22, negleda32 = mapa[key].split()

	if strand != strand2:
		zastavica = False
		if key in dupliciSet:
			for line2 in dupliciList:
				if (line2.split()[0],line2.split()[5]) == key and line2.split()[4] == strand:
					zastavica = True
					break
		if zastavica:
			continue
		brojRazlicitStrand += 1
		continue

	if ime == ime2 and velicina == velicina2 and pocetakQ == pocetakQ2 and krajQ == krajQ2 and strand == strand2 and REFime == REFime2 and pocetakR == pocetakR2 and krajR == krajR2:
		brojIstih += 1
		barNesto = True

	if test3(pocetakQ, pocetakQ2) and test3(krajQ, krajQ2) and test3(pocetakR, pocetakR2) and test3(krajR, krajR2):
		brojUnutarGranice3 += 1
		barNesto = True

	if test(pocetakQ, pocetakQ2) and test(krajQ, krajQ2) and test(pocetakR, pocetakR2) and test(krajR, krajR2):
		brojUnutarGranice25 += 1
		barNesto = True

	if test(pocetakQ, pocetakQ2) and test(krajQ, krajQ2):
		barNesto = True

	if test(pocetakR, pocetakR2) and test(krajR, krajR2):
		barNesto = True

	if test(pocetakQ, pocetakQ2) and test(pocetakR, pocetakR2):
		barNesto = True

	if test(krajQ, krajQ2) and test(krajR, krajR2):
		barNesto = True

	if test(pocetakQ, pocetakQ2):
		barNesto = True

	if test(krajQ, krajQ2):
		barNesto = True

	if test(pocetakR, pocetakR2):
		barNesto = True

	if test(krajR, krajR2):
		barNesto = True

	if barNesto:
		continue

	brojNemaNiceg += 1

if error != 0:
	print "Error", error

print "Minimap2 mapiranja:", velicinaMinimapa
print "Visestruka mapiranja:", duplici
print ""

print "Ukupno mapiranja:", velicinaMappera
print "Identicna mapiranja:", brojIstih
print "Unutar granice 3:", brojUnutarGranice3
print "Unutar granice 25:", brojUnutarGranice25
print "Ne pronade nista:", brojNemaNiceg
print ""

print "Ukupno mapiranja:", velicinaMappera
print "Razlicit lanac:", brojRazlicitStrand
print "Minimap2 nema:", brojMinimapNema
print "Broj losih mapiranja:", brojNemaPocetakNiKraj






