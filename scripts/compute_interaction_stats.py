
import csv
import statistics

filepath = "data/20251013_entero_busted_MSS_SRV_interactions.csv"

delta_srv = []
delta_mss = []
delta_combined = []
interactions = []

with open(filepath, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        try:
            ds = float(row['diff__BUSTED-SRV-E_minus_BUSTED-noSRV-E'])
            dm = float(row['diff__BUSTED-MSS-2classGA-noSRV_minus_BUSTED-noSRV-E'])
            dc = float(row['diff__BUSTED-MSS-2classGA-SRV_minus_BUSTED-noSRV-E'])
            
            i = (ds + dm) - dc
            
            delta_srv.append(ds)
            delta_mss.append(dm)
            delta_combined.append(dc)
            interactions.append(i)
        except (ValueError, KeyError):
            continue

n = len(interactions)
pos_i = [val for val in interactions if val > 0]
fraction_pos = len(pos_i) / n if n > 0 else 0

print(f"Total genes: {n}")
print(f"Mean Interaction (I): {statistics.mean(interactions):.2f}")
print(f"Median Interaction (I): {statistics.median(interactions):.2f}")
print(f"Fraction with I > 0 (Redundancy): {fraction_pos:.1%}")
print(f"Mean Delta LL SRV: {statistics.mean(delta_srv):.2f}")
print(f"Mean Delta LL MSS: {statistics.mean(delta_mss):.2f}")
print(f"Mean Delta LL Combined: {statistics.mean(delta_combined):.2f}")
