import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableColumn;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class frmDiffExpressedGeneListMaker {

    //region Variables
    private JPanel MainForm;
    private JButton btnCreateList;
    private JLabel lblChooseDataDirectory;
    private JButton btnChooseDataDirectory;
    private JTextField tfDataDirectory;
    private JTable tblDataCategories;
    private JButton btnAddCategory;

    // TODO I shouldn't be allowed to pass in null data...
    private GeneListParameters geneParameters = new GeneListParameters(null, new ArrayList<SampleCategory>());
    //endregion

    //region Methods
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
        initCategoryTable();
    }

    private void initCategoryTable() {
        DefaultTableModel model = new DefaultTableModel();
        model.addColumn("CEL File");
        TableColumn celFileColumn = new TableColumn(0);
        celFileColumn.setHeaderValue("CEL File");


        tblDataCategories.setModel(model);
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
                        geneParameters.CELFileDirectory = selectedFolder;
                        // TODO Setting the text as done here should be done automatically on property change.
                        tfDataDirectory.setText(geneParameters.CELFileDirectory.toString());
                        PopulateDataCategoryTable();
                    }
                }
            }
        });

        btnAddCategory.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                dlgAddCategory addCategoryDialog = new dlgAddCategory();

                // TODO This code is super hacky. Clean it up.
                if (addCategoryDialog.Display() == dlgAddCategory.DialogResult.OK) {
                    // TODO Adding to this model should force the table to update, rather than me doing it manually here.
                    geneParameters.SampleCategories.add(addCategoryDialog.ProvidedCategory);

                    ((DefaultTableModel)tblDataCategories.getModel()).addColumn(addCategoryDialog.ProvidedCategory.Name);
                    tblDataCategories.getColumnModel().getColumn(tblDataCategories.getColumnCount()-1).setCellEditor(new DefaultCellEditor(new JComboBox(addCategoryDialog.ProvidedCategory.ValidValues.toArray())));

                    addCategoryDialog.dispose();
                }
            }
        });

        btnCreateList.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                GeneListMaker.MakeGeneList(geneParameters);
            }
        });
    }

    private void PopulateDataCategoryTable() {
        for (String celFile : geneParameters.CELFileNames()) {
            ((DefaultTableModel) tblDataCategories.getModel()).addRow(new Object[] {celFile});
        }
    }
    //endregion
}