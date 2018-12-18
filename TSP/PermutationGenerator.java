package TSP;

import java.util.Arrays;
import java.util.Iterator;

/**
 * @author dtortola
 *			edited by Steven Lawrence
 *			
 *			A very nice implementation of a generator for Java.
 *         this implementation is based in Steinhaus Johnson Trotter algorithm
 *         and Shimon Even's improvement;
 * 
 * @see https://en.wikipedia.org/wiki/Steinhaus%E2%80%93Johnson%E2%80%93Trotter_algorithm
 *
 */
public class PermutationGenerator implements Iterator<int[]> {
	/**
	 * direction[i] = -1 if the element i has to move to the left, +1 to the right,
	 * 0 if it does not need to move
	 */
	private int[] direction;
	/**
	 * inversePermutation[i] is the position of element i in permutation; It's
	 * called inverse permutation because if p2 is the inverse permutation of p1,
	 * then p1 is the inverse permutation of p2
	 */
	private int[] inversePermutation;
	/**
	 * current permutation
	 */
	private int[] perm;
	/**
	 * specifies an element that doesn't move during the permutations and its index
	 */
	private int nailDown;
	private int nDownAtIndex;
	/**
	 * The number of elements involved in the permutation
	 */
	private int permCount;
	/**
	 * In case of a nail down indicates when to swap positions back
	 */
	boolean switchBack;
	

	public PermutationGenerator(int n, int nailDownPos) {
		// initialize the permutations, each integer will represent an element of the
		// set.
		perm = new int[n];
		for (int i = 0; i < n; i++) {
			perm[i] = i;
		}
		// the support elements
		inversePermutation = Arrays.copyOf(perm, n);
		direction = new int[n];
		Arrays.fill(direction, -1);
		direction[0] = 0;
		
		// take care of a nail downed position
		switchBack = false;
		
		if (nailDownPos > n || nailDownPos < 0) {
			nailDown = -1; // we don't nail down any of the positions
			permCount = n;
		} else {
			nailDown = nailDownPos;
			permCount = n - 1;
		}

	}

	public PermutationGenerator(int n) {
		this(n, -1);
	}

	/**
	 * Swaps the elements in array at positions i1 and i2
	 * 
	 * @param array
	 * @param i1
	 * @param i2
	 */
	private static void swap(int[] array, int i1, int i2) {
		int temp = array[i1];
		array[i1] = array[i2];
		array[i2] = temp;
	}

	/**
	 * prepares permutation to be the next one to return
	 */
	private void buildNextPermutation() {
		
		
		if (switchBack) {
			swap(inversePermutation, perm[nailDown], perm[nDownAtIndex]);
			swap(perm, nailDown, nDownAtIndex);
			
			swap(inversePermutation, perm[nailDown], perm[perm.length - 1]);
			swap(perm, nailDown, perm.length - 1);
			
			switchBack = false;
		}
		
		// find the largest element with a nonzero direction, and swaps it in
		// the indicated direction
		int index = -1;
		for (int i = 0; i < permCount; i++) {
			if (direction[perm[i]] != 0 && (index < 0 || perm[index] < perm[i])) {
				index = i;
			}
		}
		if (index < 0) {
			// there are no more permutations
			perm = null;
		} else {
			
			// element we're moving
			int chosenElement = perm[index];
			// direction we're moving
			int dir = direction[chosenElement];
			// index2 is the new position of chosenElement
			int index2 = index + dir;

			// we'll swap positions elements permutation[index] and
			// permutation[index2] in permutation, to keep inversePermutation we
			// have to swap inversePermutation's elements at index
			// permutation[index] and permutation[index2]
			/*
			 * System.out.println(index); System.out.println(index2);
			 * System.out.println(perm[index]); System.out.println(direction[index]);
			 */
			
			swap(inversePermutation, perm[index], perm[index2]);
			swap(perm, index, index2);
			

			// update directions
			if (index2 == 0 || index2 == permCount - 1 || perm[index2 + dir] > perm[index2]) {
				// direction of chosen element
				direction[chosenElement] = 0;
			}

			// all elements greater that chosenElement set its direction to +1
			// if they're before index-1 or -1 if they're after
			for (int i = chosenElement + 1; i < permCount; i++) {
				if (inversePermutation[i] > index2) {
					direction[i] = -1;
				} else {
					direction[i] = 1;
				}
			}
			
			// handle nail down
			if (nailDown != -1 && perm[nailDown] != nailDown ) {
				
				swap(inversePermutation, perm[nailDown], perm[perm.length - 1]);
				swap(perm, nailDown, perm.length - 1);
				
				nDownAtIndex = 0;
				for (int in = 0; in < perm.length; in++) {
					if (perm[in] == nailDown)
						nDownAtIndex = in;
				}
				swap(inversePermutation, perm[nailDown], perm[nDownAtIndex]);
				swap(perm, nailDown, nDownAtIndex);
				
				switchBack = true;
			}
		}
	}
	
	// stops when their is no more permutations to do
	@Override
	public boolean hasNext() {
		return perm != null;
	}

	@Override
	public int[] next() {
		int[] result = Arrays.copyOf(perm, perm.length);
		buildNextPermutation();
		return result;
	}
	
	// testing
	public static void main(String[] args) {
		int[] gotit;
		int n = 0;
		PermutationGenerator pg = new PermutationGenerator(6);
		while (pg.hasNext()) {
			gotit = pg.next();
			for (int i = 0; i < gotit.length; i++) {
				System.out.print(gotit[i] + " ");
			}
			System.out.println();

			n++;
		}
		System.out.println(n);
	}
	
	public void printArray(int[] p, String mes) {
		System.out.print(mes);
		String s = "";
		for (int b : p) {
			s += b + " ";
		}
		System.out.println(s);
	}
	
	public void printArray(int[] p) {
		printArray(p, "");
	}

}
