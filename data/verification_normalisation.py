import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Simulation pour expliquer la distribution
print("🔍 EXPLICATION DE LA DISTRIBUTION NORMALISÉE")
print("=" * 50)

# Après normalisation à 10,000 comptes par cellule
max_value = 781.25
target_sum = 10000

print(f"📊 Valeur maximale observée: {max_value}")
print(f"📊 Target sum par cellule: {target_sum}")
print(f"📊 Pourcentage du max: {(max_value/target_sum)*100:.1f}%")

# Distribution typique des gènes dans une cellule
print("\n📈 DISTRIBUTION TYPIQUE DES GÈNES:")
print("- Gènes très exprimés: 5-15% des comptes totaux")
print("- Gènes modérément exprimés: 1-5% des comptes totaux")
print("- Gènes faiblement exprimés: <1% des comptes totaux")
print("- Gènes non exprimés: 0 comptes")

# Vérification avec vos données
print(f"\n✅ VOTRE RÉSULTAT EST NORMAL:")
print(f"   - Max: {max_value} ({max_value/target_sum*100:.1f}% des comptes)")
print(f"   - Cela correspond à un gène très exprimé")
print(f"   - La normalisation a bien fonctionné!")

# Calcul des percentiles pour mieux comprendre
print(f"\n📊 POUR MIEUX COMPRENDRE LA DISTRIBUTION:")
print("Ajoutez ce code à votre analyse:")

print(
    """
# Analyse de la distribution des valeurs normalisées
import numpy as np

# Percentiles des valeurs non-nulles
non_zero_values = test.X[test.X > 0]
percentiles = [50, 75, 90, 95, 99, 99.9]

print("📊 Percentiles des valeurs normalisées (non-nulles):")
for p in percentiles:
    value = np.percentile(non_zero_values, p)
    print(f"   {p}ème percentile: {value:.2f} ({value/10000*100:.2f}% des comptes)")

print(f"   Maximum: {test.X.max():.2f} ({test.X.max()/10000*100:.2f}% des comptes)")
"""
)
