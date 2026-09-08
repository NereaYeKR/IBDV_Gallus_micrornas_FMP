import sys

def rename_clusters_mature():
    # Nombres de archivos según tu estructura actual
    clusters_file = 'mature-clusters.txt'
    counts_file = 'mature-counts.txt'
    output_file = 'mature-counts-gga.txt'

    # 1. Crear el diccionario de mapeo
    # Un Cluster-ID puede tener varios gga-miR asociados
    mapping = {}
    
    print(f"Leyendo mapeo desde {clusters_file}...")
    try:
        with open(clusters_file, 'r') as f:
            for line in f:
                line = line.strip()
                # Saltamos líneas vacías o con formato extraño
                if not line or line.startswith('['): 
                    continue
                
                # Usamos split() sin argumentos para que detecte cualquier 
                # cantidad de espacios o tabuladores
                parts = line.split()
                if len(parts) >= 2:
                    mir_id = parts[0]
                    cluster_id = parts[1]
                    
                    if cluster_id not in mapping:
                        mapping[cluster_id] = []
                    mapping[cluster_id].append(mir_id)
    except FileNotFoundError:
        print(f"Error: No se encuentra el archivo {clusters_file}")
        return

    # 2. Procesar el archivo de counts y duplicar filas
    print(f"Procesando counts desde {counts_file}...")
    try:
        with open(counts_file, 'r') as f_in, open(output_file, 'w') as f_out:
            # Leer la cabecera
            header = f_in.readline()
            f_out.write(header)
            
            count_processed = 0
            count_expanded = 0
            
            for line in f_in:
                if not line.strip(): continue
                
                parts = line.strip().split('\t')
                cluster_name = parts[0]
                
                # Si el cluster tiene mapeo a uno o más miRNAs
                if cluster_name in mapping:
                    for mir_id in mapping[cluster_name]:
                        # Creamos la nueva línea: miR_ID + el resto de las columnas de counts
                        new_line = mir_id + '\t' + '\t'.join(parts[1:]) + '\n'
                        f_out.write(new_line)
                        count_expanded += 1
                else:
                    # Si no hay mapeo, mantenemos la línea original por seguridad
                    f_out.write(line)
                
                count_processed += 1
        
        print("-" * 30)
        print(f"¡Éxito!")
        print(f"Clusters procesados: {count_processed}")
        print(f"Filas totales en el nuevo archivo: {count_expanded}")
        print(f"Archivo generado: '{output_file}'")
                
    except FileNotFoundError:
        print(f"Error: No se encuentra el archivo {counts_file}")

if __name__ == "__main__":
    rename_clusters_mature()