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
brojLosaVelicina = 0

for line in content:
	key = line.split()[0]

	if int(line.split()[1]) < (2 * k_value):
		brojLosaVelicina += 1

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
brojUnutarGraniceQuery = 0
brojUnutarGraniceReference = 0
brojRazlicitStrand = 0
brojUnutarGranicePocetakQuery = 0
brojUnutarGraniceKrajQuery = 0
brojUnutarGranicePocetakReference = 0
brojUnutarGraniceKrajReference = 0
brojUnutarGranicePocetakOba = 0
brojUnutarGraniceKrajOba = 0
brojNemaNiceg = 0
error = 0

for line in content:
	key, velicina, pocetakQ, krajQ, strand, REFime, REFvelicina, pocetakR, krajR, negleda1, negleda2, negleda3 = line.split()
	

	barNesto = False

	if int(krajQ) > int(velicina):
		error += 1
		continue

	if key not in mapa:
		brojMinimapNema += 1
		continue

	key2, velicina2, pocetakQ2, krajQ2, strand2, REFime2, REFvelicina2, pocetakR2, krajR2, negleda12, negleda22, negleda32 = mapa[key].split()

	if strand != strand2:
		zastavica = False
		if key in dupliciSet:
			for line2 in dupliciList:
				if line2.split()[0] == key and line2.split()[4] == strand:
					zastavica = True
					break
		if zastavica:
			continue
		brojRazlicitStrand += 1
		continue

	if key == key2 and velicina == velicina2 and pocetakQ == pocetakQ2 and krajQ == krajQ2 and strand == strand2 and REFime == REFime2 and pocetakR == pocetakR2 and krajR == krajR2:
		brojIstih += 1
		barNesto = True

	if test3(pocetakQ, pocetakQ2) and test3(krajQ, krajQ2) and test3(pocetakR, pocetakR2) and test3(krajR, krajR2):
		brojUnutarGranice3 += 1
		barNesto = True

	if test(pocetakQ, pocetakQ2) and test(krajQ, krajQ2) and test(pocetakR, pocetakR2) and test(krajR, krajR2):
		brojUnutarGranice25 += 1
		barNesto = True

	if test(pocetakQ, pocetakQ2) and test(krajQ, krajQ2):
		brojUnutarGraniceQuery += 1
		barNesto = True

	if test(pocetakR, pocetakR2) and test(krajR, krajR2):
		brojUnutarGraniceReference +=1
		barNesto = True

	if test(pocetakQ, pocetakQ2) and test(pocetakR, pocetakR2):
		brojUnutarGranicePocetakOba +=1
		barNesto = True

	if test(krajQ, krajQ2) and test(krajR, krajR2):
		brojUnutarGraniceKrajOba +=1
		barNesto = True

	if test(pocetakQ, pocetakQ2):
		brojUnutarGranicePocetakQuery +=1
		barNesto = True

	if test(krajQ, krajQ2):
		brojUnutarGraniceKrajQuery +=1
		barNesto = True

	if test(pocetakR, pocetakR2):
		brojUnutarGranicePocetakReference +=1
		barNesto = True

	if test(krajR, krajR2):
		brojUnutarGraniceKrajReference +=1
		barNesto = True

	if barNesto:
		continue

	brojNemaNiceg += 1


	#missPP = int(pocetakQ2) - int(pocetakQ)
	#missPK = int(krajQ2) - int(krajQ)
	#missKP = int(pocetakR2) - int(pocetakR)
	#missKK = int(krajR2) - int(krajR)
	
	#print key,"\t", velicina, "\t", missPP, "\t", missPK, "\t", strand, "\t", missKP, "\t", missKK
	#print line
	#print mapa[key]
	#print "\n"

if error != 0:
	print "Error", error

print "Broj dobiven minimapom", velicinaMinimapa
print "Duplici u minimapu", duplici
print "Broj sekvenci s duljinom manjom od 2K", brojLosaVelicina
print "Los pocetak i kraj", brojNemaPocetakNiKraj
print ""
print "Broj dobivem mapperom", velicinaMappera
print "Broj zapisa kojih minimap nema", brojMinimapNema
print ""
print "Broj razlicit strand", brojRazlicitStrand
print "Broj istih", brojIstih
print "Broj unutar granice 3:", brojUnutarGranice3
print "Broj unutar granice 25:", brojUnutarGranice25
print "Broj unutar granice (samo query)", brojUnutarGraniceQuery
print "Broj unutar granice (samo reference)", brojUnutarGraniceReference
print "Broj unutar granice (pocetak oba)", brojUnutarGranicePocetakOba
print "Broj unutar granice (kraj oba)", brojUnutarGraniceKrajOba
print "Broj unutar granice (query pocetak)", brojUnutarGranicePocetakQuery
print "Broj unutar granice (query kraj)", brojUnutarGraniceKrajQuery
print "Broj unutar granice (reference pocetak)", brojUnutarGranicePocetakReference
print "Broj unutar granice (reference kraj)", brojUnutarGraniceKrajReference
print "Broj nije nadeno nista", brojNemaNiceg







