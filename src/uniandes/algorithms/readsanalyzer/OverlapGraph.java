package uniandes.algorithms.readsanalyzer;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;
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
		if (readCounts.containsKey(sequence)) {
			readCounts.compute(sequence, (key, i) -> i + 1);
		} else {
			readCounts.put(sequence, 1);
			ArrayList<ReadOverlap> sequenceOverlaps = new ArrayList<>();
			for (String sequence2 : getDistinctSequences()) {
				int overlap = getOverlapLength(sequence, sequence2);
				if (overlap >= minOverlap)
					sequenceOverlaps.add(new ReadOverlap(sequence, sequence2, overlap));
			}
			overlaps.put(sequence, sequenceOverlaps);

			for (Entry<String, ArrayList<ReadOverlap>> entry : overlaps.entrySet()) {
				int overlap = getOverlapLength(entry.getKey(), sequence);
				if (overlap >= minOverlap)
					entry.getValue().add(new ReadOverlap(entry.getKey(), sequence, overlap));
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
		ext: for (int length = Math.min(sequence2.length(), sequence1.length()); length > 0; length--) {
			
			for (int i = sequence1.length() - length, j = 0; i < sequence1.length(); i++, j++) {
				sequence1.charAt(i);
				sequence2.charAt(j);
				if (sequence1.charAt(i) != sequence2.charAt(j))
					continue ext;
			}
				
				
			return length;
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
		int[] abundancias = new int[readCounts.values().stream().max(Integer::compare).get() + 1];
		 readCounts.values().stream().forEach(conteo -> abundancias[conteo]++);
		return abundancias;
	}

	/**
	 * Calculates the distribution of number of successors
	 * 
	 * @return int [] array where the indexes are number of successors and the
	 *         values are the number of sequences having as many successors as the
	 *         corresponding array index.
	 */
	public int[] calculateOverlapDistribution() {
		int[] sucesores = new int[ overlaps.values().stream().map((array) -> array.size()).max(Integer::compare).get() + 1];
		 overlaps.values().stream().map((array) -> array.size()).forEach((numero_sucesores) -> sucesores[numero_sucesores]++);
		return sucesores;
	}

	/**
	 * Predicts the leftmost sequence of the final assembly for this overlap graph
	 * 
	 * @return String Source sequence for the layout path that will be the left most
	 *         subsequence in the assembly
	 */
	public String getSourceSequence() {
		return overlaps.values().stream()
				.flatMap(Collection::stream) // Aggregate sub lists to one stream
				.collect(Collectors.groupingBy(ReadOverlap::getDestSequence))// group by destSequence 
				.entrySet().stream()
				.min((a, b) -> a.getValue().size() - b.getValue().size()) // Minimum count
				.get().getKey(); // Get String
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
		String act = getSourceSequence();
		while (true) {
			visitedSequences.add(act);
			Optional<ReadOverlap> next = overlaps.get(act).stream()
					.filter(a -> !visitedSequences.contains(a.getDestSequence()))
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
		assembly.append(layout.get(0).getSourceSequence());
		for (ReadOverlap r_layout : layout)
			assembly.append(r_layout.getDestSequence().substring(r_layout.getOverlap()));
		return assembly.toString();
	}

}
