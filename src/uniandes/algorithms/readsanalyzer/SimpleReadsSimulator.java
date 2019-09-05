package uniandes.algorithms.readsanalyzer;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;
import ngsep.sequences.io.FastaSequencesHandler;

/**
 * Simple script that simulates error free reads from a text in fasta format
 * 
 * @author Jorge Duitama
 *
 */
public class SimpleReadsSimulator {
	/**
	 * Main class that executes the program
	 * 
	 * @param args Array of arguments: args[0]: Source sequence in fasta format. If
	 *             many sequences are present, it only takes the first sequence
	 *             args[1]: Length of the reads to simulate args[2]: Number of reads
	 *             to simulate args[3]: Path to the output file
	 * @throws Exception If the fasta file can not be loaded
	 */
	public static void main(String[] args) throws Exception {
		String filename = args[0];
		int readLength = Integer.parseInt(args[1]);
		int numReads = Integer.parseInt(args[2]);
		String outFile = args[3];
		FastaSequencesHandler handler = new FastaSequencesHandler();
		handler.setSequenceType(StringBuilder.class);
		QualifiedSequenceList sequences = handler.loadSequences(filename);
		if (sequences.size() == 0)
			throw new Exception("No sequences found in file: " + filename);
		QualifiedSequence seq = sequences.get(0);
		String sequence = seq.getCharacters().toString();
		int seqLength = sequence.length();
		System.out.println("Length of the sequence to simulate reads: " + seqLength);
		double averageRD = ((double) numReads * readLength) / seqLength;
		System.out.println("Expected average RD: " + averageRD);
		char[] fixedQS = new char[readLength];
		Arrays.fill(fixedQS, '5');
		String fixedQSStr = new String(fixedQS);
		Random random = new Random();
		double tasaCambios = Double.valueOf(args[4]);
		double tasaIndels = Double.valueOf(args[5]);

		try (PrintStream out = new PrintStream(outFile)) {
			for (int i = 0; i < numReads; i++) {
				out.println(">" +seq.getName()+" " + seq.getComments());
				int pos = random.nextInt(seqLength - readLength + 1);
				String sub = sequence.substring(pos, pos + readLength);
				
				sub = agregarIndels(random, tasaIndels, sub);
				sub = agregarErrores(random, tasaCambios, sub);
				out.println(sub);
				/**out.println("+");
				out.println(fixedQSStr);*/
			}
		}
	}

	private static String agregarIndels(Random random, double tasaIndels, String sub) {
		TreeSet<Integer> positions = new TreeSet<Integer>();
		for (int j = 0; j < sub.length() * tasaIndels; j++) {
			int p = random.nextInt(sub.length());
			while (positions.contains(p))
				p = random.nextInt(sub.length());
			positions.add(p);
		}
		char[] ans = new char[sub.length() - positions.size()];
		int l = 0, j = 0;
		for (int skip : positions) {
			while (l < skip)
				ans[l++] = sub.charAt(j++);
			j++;
		}
		while (j < sub.length())
			ans[l++] = sub.charAt(j++);

		return new String(ans);
	}

	private static String agregarErrores(Random random, double tasaCambios, String sub) {
		char[] letras = new char[] { 'A', 'C', 'T', 'G' };
		char[] array = sub.toCharArray();
		Set<Integer> changedIndex = new HashSet<>();
		for (int j = 0; j < array.length * tasaCambios; j++) {
			int p = random.nextInt(array.length);
			while (changedIndex.contains(p))
				p = random.nextInt(array.length);
			changedIndex.add(p);

			char let = letras[random.nextInt(letras.length)];
			while (let == array[p])
				let = letras[random.nextInt(letras.length)];
		}
		return new String(array);
	}
}
