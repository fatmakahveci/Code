import java.util.ArrayList;
import java.util.List;

public class MergeTwoSortedLists {
    public static void main(String[] args) {
        List<Integer> firstList=new ArrayList<Integer>();
        List<Integer> secondList=new ArrayList<Integer>();
        List<Integer> mergedList=new ArrayList<Integer>();
        firstList.add(1);
        firstList.add(3);
        firstList.add(4);
        firstList.add(5);
        firstList.add(9);
        firstList.add(10);

        secondList.add(2);
        secondList.add(4);
        secondList.add(5);
        secondList.add(6);
        secondList.add(7);
        secondList.add(11);

        int firstListIndex=0;
        int secondListIndex=0;

        while(firstListIndex < firstList.size() && secondListIndex < secondList.size()) {
            if(firstList.get(firstListIndex) <= secondList.get(secondListIndex)) {
                mergedList.add(firstList.get(firstListIndex));
                firstListIndex++;
            } else {
                mergedList.add(secondList.get(secondListIndex));
                secondListIndex++;
            }
        }
        while(firstListIndex<firstList.size()) {
            mergedList.add(firstList.get(firstListIndex));
            firstListIndex++;
        }

        while(secondListIndex<secondList.size()) {
            mergedList.add(secondList.get(secondListIndex));
            secondListIndex++;
        }
        System.out.println(mergedList.toString());
    }
}
