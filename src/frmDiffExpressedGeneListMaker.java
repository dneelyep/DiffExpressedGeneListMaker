import javax.swing.*;

public class frmDiffExpressedGeneListMaker {

    private JPanel MainForm;

    public static void main(String[] args) {
        JFrame frame = new JFrame("frmDiffExpressedGeneListMaker");
        frame.setContentPane(new frmDiffExpressedGeneListMaker().MainForm);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.pack();
        frame.setVisible(true);
    }
}