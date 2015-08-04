import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

public class frmDiffExpressedGeneListMaker {

    private JPanel MainForm;
    private JButton btnCreateList;
    private JLabel lblChooseDataDirectory;
    private JButton btnChooseDataDirectory;
    private JTextField tfDataDirectory;

    private GeneListParameters geneParameters;


    public static void main(String[] args) {
        JFrame frame = new JFrame("DEGLMaker (Differentially-expressed gene list maker)");
        frame.setContentPane(new frmDiffExpressedGeneListMaker().MainForm);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.pack();

        // Center the form on-screen.
        frame.setLocationRelativeTo(null);

        frame.setVisible(true);
    }

    public frmDiffExpressedGeneListMaker() {
        initActionListeners();
    }

    private void initActionListeners() {
        btnChooseDataDirectory.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                JFileChooser fc = new JFileChooser();
                fc.setDialogTitle("Open data directory");
                fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                int returnVal = fc.showOpenDialog(MainForm);

                if (returnVal == JFileChooser.APPROVE_OPTION) {
                    File selectedFolder = fc.getSelectedFile();
                    if (selectedFolder != null) {
                        geneParameters = new GeneListParameters(selectedFolder);
                        // TODO Setting the text as done here should be done automatically on property change.
                        tfDataDirectory.setText(geneParameters.CELFileDirectory.toString());
                    }
                }
            }
        });
    }
}