def generate_pdb(center, size):
    # Extract center and size
    cx, cy, cz = center
    sx, sy, sz = size

    # Calculate the half sizes
    hx, hy, hz = sx / 2, sy / 2, sz / 2

    # Define the vertices based on center and half sizes
    vertices = [
        (cx - hx, cy - hy, cz - hz),
        (cx - hx, cy - hy, cz + hz),
        (cx - hx, cy + hy, cz - hz),
        (cx + hx, cy - hy, cz - hz),
        (cx + hx, cy - hy, cz + hz),
        (cx + hx, cy + hy, cz - hz),
        (cx - hx, cy + hy, cz + hz),
        (cx + hx, cy + hy, cz + hz),
    ]

    # Define the center point
    center_point = (cx, cy, cz)

    # Start writing the PDB content
    pdb_content = []

    # Write the vertices
    for i, vertex in enumerate(vertices, start=1):
        pdb_content.append(
            f"HETATM{str(i).rjust(5)} PS{i}  PSD P   1    {vertex[0]:8.3f}{vertex[1]:8.3f}{vertex[2]:8.3f}  0.00  0.00          PSDOPS"
        )

    # Write the center point
    pdb_content.append(
        f"HETATM{str(len(vertices) + 1).rjust(5)} PS{len(vertices) + 1}  PSD P   1    {center_point[0]:8.3f}{center_point[1]:8.3f}{center_point[2]:8.3f}  0.00  0.00          PSDOPS"
    )

    # Write the CONECT records
    conect_records = [
        "CONECT    1    2    4    3",
        "CONECT    2    1    5    7",
        "CONECT    3    1    6    7",
        "CONECT    4    1    5    6",
        "CONECT    5    2    4    8",
        "CONECT    6    3    4    8",
        "CONECT    7    2    3    8",
        "CONECT    8    5    6    7",
    ]

    pdb_content.extend(conect_records)
    pdb_content.append("END")

    # Join the PDB content into a single string
    pdb_string = "\n".join(pdb_content)

    return pdb_string


# Example usage
center = (0, 0, 0)
size = (8, 4, 6)
pdb_output = generate_pdb(center, size)
print(pdb_output)

# To save the output to a file:
with open("grid_box.pdb", "w") as f:
    f.write(pdb_output)
