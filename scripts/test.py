print("\n===== GLOBAL FLOW CHECK =====")
print("Total target bait capacity:", sum(target_bait))
print("Total target prey capacity:", sum(target_prey))
print("Theoretical max flow:", min(sum(target_bait), sum(target_prey)))
print("Actual flow:", flow_value)
print("Flow / theoretical max:",
      flow_value / max(1, min(sum(target_bait), sum(target_prey))))

neg_full_g = generate_graph(negative_bait_prey_df, node_map)
neg_bait_deg, neg_prey_deg = get_degree(neg_full_g)

print("\n===== NEGATIVE GRAPH STATS =====")
print("Total negative edges:", negative_bait_prey_df.shape[0])
print("Unique negative baits:", negative_bait_prey_df['bait'].nunique())
print("Unique negative prey:", negative_bait_prey_df['prey'].nunique())
print("Mean negative bait degree:", np.mean(neg_bait_deg))
print("Mean negative prey degree:", np.mean(neg_prey_deg))

print("\n===== BAIT BOTTLENECK CHECK =====")

bait_limited = 0
for i in range(len(target_bait)):
    if target_bait[i] > neg_bait_deg[i]:
        bait_limited += 1

print("Baits where target > available negative degree:", bait_limited)
print("Fraction bait-limited:",
      bait_limited / len(target_bait))


print("\n===== PREY BOTTLENECK CHECK =====")

prey_limited = 0
for i in range(len(target_prey)):
    if target_prey[i] > neg_prey_deg[i]:
        prey_limited += 1

print("Prey where target > available negative degree:", prey_limited)
print("Fraction prey-limited:",
      prey_limited / len(target_prey))

max_possible_from_baits = sum(
    min(target_bait[i], neg_bait_deg[i])
    for i in range(len(target_bait))
)

max_possible_from_prey = sum(
    min(target_prey[i], neg_prey_deg[i])
    for i in range(len(target_prey))
)

print("\n===== STRUCTURAL UPPER BOUNDS =====")
print("Upper bound from bait side:", max_possible_from_baits)
print("Upper bound from prey side:", max_possible_from_prey)
print("Actual flow:", flow_value)