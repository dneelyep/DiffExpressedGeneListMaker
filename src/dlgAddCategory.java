import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import java.awt.event.*;
import java.util.HashSet;
import java.util.Set;

public class dlgAddCategory extends JDialog {

    private JPanel contentPane;
    private JButton buttonOK;
    private JButton buttonCancel;
    private JTextField txtCategoryName;
    private JList lstCategoryValues;
    private JButton btnAdd;
    private JButton btnRemove;
    private JLabel lblCategoryName;
    private JTable tblCategoryValues;

    public dlgAddCategory() {
        setContentPane(contentPane);
        setModal(true);
        getRootPane().setDefaultButton(buttonOK);
        this.setTitle("Add experimental category");

        initCategoryModel();
        initActionListeners();


// call onCancel() when cross is clicked
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                onCancel();
            }
        });

// call onCancel() on ESCAPE
        contentPane.registerKeyboardAction(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                onCancel();
            }
        }, KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
    }

    // TODO Better name than not sure here.
    public enum DialogResult {
        NOTSURE,
        OK,
        CANCEL
    }

    public DialogResult Result = DialogResult.NOTSURE;

    public SampleCategory ProvidedCategory;

    private void initCategoryModel() {
        DefaultTableModel model = new DefaultTableModel();
        model.addColumn("Value");
        tblCategoryValues.setModel(model);
    }

    private void initActionListeners() {
        buttonOK.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                onOK();
            }
        });

        buttonCancel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                onCancel();
            }
        });

        btnAdd.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // TODO Ideally the row would be highlighted and the user could start editing immediately after pressing "Add".
                ((DefaultTableModel)tblCategoryValues.getModel()).addRow(new Object[] {"TODO Value here"});
            }
        });

        // TODO Fix form scaling positioning.
        // TODO Disable the OK button until needed info is filled in.
        btnRemove.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                int[] categoryValues = tblCategoryValues.getSelectedRows();

                for (int categoryValue : categoryValues) {
                    ((DefaultTableModel)tblCategoryValues.getModel()).removeRow(categoryValue);
                }
            }
        });
    }

    private void onOK() {
        // TODO Just write a function for converting from the UI to an object as I'm doing here.
        HashSet<String> validValues = new HashSet<>();
        for (int value = 0; value < tblCategoryValues.getRowCount(); value++) {
            validValues.add((String) tblCategoryValues.getValueAt(value, 0));
        }

        ProvidedCategory = new SampleCategory(txtCategoryName.getText(),
                                              validValues);

        Result = DialogResult.OK;
        dispose();
    }

    private void onCancel() {
        Result = DialogResult.CANCEL;
        dispose();
    }

    public DialogResult Display() {
        pack();

        // TODO Make an extension method or similar to center a form.
        // TODO This is not working.
        // Center the dialog on-screen.
        setLocationRelativeTo(null);

        setVisible(true);

        return Result;
    }
}