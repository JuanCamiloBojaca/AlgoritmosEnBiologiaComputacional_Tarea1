package uniandes.algorithms.readsanalyzer;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import ngsep.sequences.RawRead;

/**
 * Stores abundances information on a list of subsequences of a fixed length k
 * (k-mers)
 * 
 * @author Jorge Duitama
 */
public class KmersTable implements RawReadProcessor {
	private int kmerSize;
	private Map<String, Integer> conteoKmer;

	/**
	 * Creates a new table with the given k-mer size
	 * 
	 * @param kmerSize length of k-mers stored in this table
	 */
	public KmersTable(int kmerSize) {
		// TODO: Implementar metodo
		this.kmerSize = kmerSize;
		conteoKmer = new HashMap<String, Integer>();
	}

	/**
	 * Identifies k-mers in the given read
	 * 
	 * @param read object to extract new k-mers
	 */
	public void processRead(RawRead read) {
		String sequence = read.getSequenceString();
		// TODO Implementar metodo. Calcular todos los k-mers del tamanho dado en la
		// constructora y actualizar la abundancia de cada k-mer
		for (int i = 0; i <= sequence.length() - kmerSize; i++) {
			String kmer = sequence.substring(i, i + kmerSize);
			if (conteoKmer.containsKey(kmer))
				conteoKmer.compute(kmer, (key, j) -> j + 1);
			else
				conteoKmer.put(kmer, 1);
		}
	}

	/**
	 * List with the different k-mers found up to this point
	 * 
	 * @return Set<String> set of k-mers
	 */
	public Set<String> getDistinctKmers() {
		// TODO Implementar metodo
		return conteoKmer.keySet();
	}

	/**
	 * Calculates the current abundance of the given k-mer
	 * 
	 * @param kmer sequence of length k
	 * @return int times that the given k-mer have been extracted from given reads
	 */
	public int getAbundance(String kmer) {
		// TODO Implementar metodo
		return conteoKmer.get(kmer);
	}

	/**
	 * Calculates the distribution of abundances
	 * 
	 * @return int [] array where the indexes are abundances and the values are the
	 *         number of k-mers observed as many times as the corresponding array
	 *         index. Position zero should be equal to zero
	 */
	public int[] calculateAbundancesDistribution() {
		// TODO Implementar metodo
		int max_length = conteoKmer.values().stream().max((a, b) -> a - b).get();
		int[] ans = new int[max_length + 1];
		conteoKmer.values().stream().forEach((conteo) -> ans[conteo]++);
		return ans;
	}
}
