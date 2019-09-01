package uniandes.algorithms.readsanalyzer;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.Set;
import ngsep.sequences.RawRead;

/**
 * Represents an overlap graph for a set of reads taken from a sequence to
 * assemble
 * 
 * @author Jorge Duitama
 *
 */
public class OverlapGraph implements RawReadProcessor {

	private int minOverlap;
	private Map<String, Integer> readCounts = new HashMap<>();
	private Map<String, ArrayList<ReadOverlap>> overlaps = new HashMap<>();

	/**
	 * Creates a new overlap graph with the given minimum overlap
	 * 
	 * @param minOverlap Minimum overlap
	 */
	public OverlapGraph(int minOverlap) {
		this.minOverlap = minOverlap;
	}

	/**
	 * Adds a new read to the overlap graph
	 * 
	 * @param read object with the new read
	 */
	public void processRead(RawRead read) {
		String sequence = read.getSequenceString();
		// TODO: Paso 1. Agregar la secuencia al mapa de conteos si no existe.
		// Si ya existe, solo se le suma 1 a su conteo correspondiente y no se deben
		// ejecutar los pasos 2 y 3
		if (readCounts.containsKey(sequence)) {
			readCounts.compute(sequence, (key, i) -> i + 1);
		} else {
			readCounts.put(sequence, 1);

			// TODO: Paso 2. Actualizar el mapa de sobrelapes con los sobrelapes en los que
			// la secuencia nueva sea predecesora de una secuencia existente
			// 2.1 Crear un ArrayList para guardar las secuencias que tengan como prefijo un
			// sufijo de la nueva secuencia
			ArrayList<ReadOverlap> lista = new ArrayList<>();
			// 2.2 Recorrer las secuencias existentes para llenar este ArrayList creando los
			// nuevos sobrelapes que se encuentren.
			for (String seq : getDistinctSequences()) {
				int sobrelape = getOverlapLength(sequence, seq);
				if (sobrelape >= minOverlap)
					lista.add(new ReadOverlap(sequence, seq, sobrelape));
			}
			// 2.3 Después del recorrido para llenar la lista, agregar la nueva secuencia
			// con su lista de sucesores al mapa de sobrelapes
			overlaps.put(sequence, lista);

			// TODO: Paso 3. Actualizar el mapa de sobrelapes con los sobrelapes en los que
			// la secuencia nueva sea sucesora de una secuencia existente
			// Recorrer el mapa de sobrelapes. Para cada secuencia existente que tenga como
			// sufijo un prefijo de la nueva secuencia
			// se agrega un nuevo sobrelape a la lista de sobrelapes de la secuencia
			// existente
			for (Entry<String, ArrayList<ReadOverlap>> a : overlaps.entrySet()) {
				int sobrelape = getOverlapLength(a.getKey(), sequence);
				if (sobrelape >= minOverlap)
					a.getValue().add(new ReadOverlap(a.getKey(), sequence, sobrelape));
			}

		}
	}

	/**
	 * Returns the length of the maximum overlap between a suffix of sequence 1 and
	 * a prefix of sequence 2
	 * 
	 * @param sequence1 Sequence to evaluate suffixes
	 * @param sequence2 Sequence to evaluate prefixes
	 * @return int Maximum overlap between a prefix of sequence2 and a suffix of
	 *         sequence 1
	 */
	private int getOverlapLength(String sequence1, String sequence2) {
		for (int i = Math.min(sequence2.length(), sequence1.length()); i > 0; i--) {
			String suffijo = sequence1.substring(sequence1.length() - i);
			String prefijo = sequence2.substring(0, suffijo.length());

			if (suffijo.equals(prefijo))
				return suffijo.length();
		}
		return 0;
	}

	/**
	 * Returns a set of the sequences that have been added to this graph
	 * 
	 * @return Set<String> of the different sequences
	 */
	public Set<String> getDistinctSequences() {
		return readCounts.keySet();
	}

	/**
	 * Calculates the abundance of the given sequence
	 * 
	 * @param sequence to search
	 * @return int Times that the given sequence has been added to this graph
	 */
	public int getSequenceAbundance(String sequence) {
		// TODO: Implementar metodo
		return readCounts.get(sequence);
	}

	/**
	 * Calculates the distribution of abundances
	 * 
	 * @return int [] array where the indexes are abundances and the values are the
	 *         number of sequences observed as many times as the corresponding array
	 *         index. Position zero should be equal to zero
	 */
	public int[] calculateAbundancesDistribution() {
		// TODO: Implementar metodo
		int max = readCounts.values().stream().max((Integer i, Integer j) -> i - j).get();
		int[] ans = new int[max + 1];
		for (int conteo : readCounts.values())
			ans[conteo]++;
		return ans;
	}

	/**
	 * Calculates the distribution of number of successors
	 * 
	 * @return int [] array where the indexes are number of successors and the
	 *         values are the number of sequences having as many successors as the
	 *         corresponding array index.
	 */
	public int[] calculateOverlapDistribution() {
		// TODO: Implementar metodo
		int max_length = overlaps.values().stream().max((a, b) -> a.size() - b.size()).get().size();
		int[] sucesores = new int[max_length + 1];
		overlaps.values().stream().mapToInt((array) -> array.size())
				.forEach((numero_sucesores) -> sucesores[numero_sucesores]++);
		return sucesores;
	}

	/**
	 * Predicts the leftmost sequence of the final assembly for this overlap graph
	 * 
	 * @return String Source sequence for the layout path that will be the left most
	 *         subsequence in the assembly
	 */
	public String getSourceSequence() {
		// TODO Implementar metodo recorriendo las secuencias existentes y buscando una
		// secuencia que no tenga predecesores
		Map<String, Integer> numeroPredecesores = new HashMap<String, Integer>();
		for (ArrayList<ReadOverlap> list : overlaps.values()) {
			for (ReadOverlap sobrelape : list) {
				String act = sobrelape.getDestSequence();
				if (numeroPredecesores.containsKey(act)) {
					numeroPredecesores.compute(act, (key, i) -> i + 1);
				} else
					numeroPredecesores.put(act, 1);
			}
		}
		return numeroPredecesores.entrySet().stream().min((a, b) -> a.getValue() - b.getValue()).get().getKey();
	}

	/**
	 * Calculates a layout path for this overlap graph
	 * 
	 * @return ArrayList<ReadOverlap> List of adjacent overlaps. The destination
	 *         sequence of the overlap in position i must be the source sequence of
	 *         the overlap in position i+1.
	 */
	public ArrayList<ReadOverlap> getLayoutPath() {
		ArrayList<ReadOverlap> layout = new ArrayList<>();
		HashSet<String> visitedSequences = new HashSet<>();
		// TODO Implementar metodo. Comenzar por la secuencia fuente que calcula el
		// método anterior

		// Luego, hacer un ciclo en el que en cada paso se busca la secuencia no
		// visitada que tenga mayor sobrelape con la secuencia actual.
		// Agregar el sobrelape a la lista de respuesta y la secuencia destino al
		// conjunto de secuencias visitadas. Parar cuando no se encuentre una secuencia
		// nueva

		String act = getSourceSequence();

		while (true) {
			visitedSequences.add(act);
			Optional<ReadOverlap> next = overlaps.get(act).stream()
					.filter((a) -> !visitedSequences.contains(a.getDestSequence()))
					.max((a, b) -> a.getOverlap() - b.getOverlap());
			if (next.isPresent()) {
				act = next.get().getDestSequence();
				layout.add(next.get());
			} else
				break;
		}

		return layout;
	}

	/**
	 * Predicts an assembly consistent with this overlap graph
	 * 
	 * @return String assembly explaining the reads and the overlaps in this graph
	 */
	public String getAssembly() {
		ArrayList<ReadOverlap> layout = getLayoutPath();
		StringBuilder assembly = new StringBuilder();
		// TODO Recorrer el layout y ensamblar la secuencia agregando al objeto assembly
		// las bases adicionales que aporta la región de cada secuencia destino que
		// está a la derecha del sobrelape
		assembly.append(layout.get(0).getSourceSequence());
		for (ReadOverlap r_layout : layout)
			assembly.append(r_layout.getDestSequence().substring(r_layout.getOverlap()));

		return assembly.toString();
	}

}
