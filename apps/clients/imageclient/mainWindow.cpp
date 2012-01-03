#include "mainWindow.h"

#include <QString>
#include <QFileDialog>
#include <QMessageBox>

#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "GadgetSocketSender.h" //We need this for now for the GadgetAcquisitionWriter
#include "GadgetMRIHeaders.h"
#include "GadgetContainerMessage.h"
#include "hoNDArray.h"
#include "../mriclient/ImageWriter.h"
#include "FileInfo.h"

const char *IMAGE_CLIENT_XML_PARAMETER_SCRIPT =
"<?xml version=\"1.0\" ?> \
   <gadgetron> \
     <encoding> \
       <matrix_size> \
          <value>192</value> \
          <value>144</value> \
       </matrix_size> \
     </encoding> \
   </gadgetron>";

//
//

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent)
{
	setupUi(this);
	retranslateUi(this);

	num_repetitions_ = 1;
	hostname_ = 0x0;
	port_no_ = 0x0;
}

void 
MainWindow::addFiles()
{
	QFileDialog fd;
	fd.setNameFilter(tr("Images (*.*)"));
	fd.setFileMode(QFileDialog::ExistingFiles);
	fd.setViewMode(QFileDialog::Detail);

	QStringList fileNames;
	if (fd.exec())
		fileNames = fd.selectedFiles();
	else return;

	listWidgetImages->addItems(fileNames);
	pushButtonGo->setEnabled(true);
}

void 
MainWindow::clearFiles()
{
	listWidgetImages->clear();
	pushButtonGo->setEnabled(false);
}

void 
MainWindow::loadConfig()
{
	QString basepath(getenv("GADGETRON_HOME"));
	if( basepath==QString("") )
		basepath = QString(".");
	else
		basepath += QString("/config");

	QFileDialog fd;
	fd.setNameFilter(tr("configs (*.xml)"));
	fd.setDirectory(QDir(basepath));
	fd.setFileMode(QFileDialog::ExistingFile);
	fd.setViewMode(QFileDialog::Detail);

	QStringList fileNames;
	if (fd.exec())
		fileNames = fd.selectedFiles();
	else return;

	QFile file(fileNames[0]);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
		return;
	QByteArray document = file.readAll();

	plainTextEditConfig->setPlainText(document);
}

void 
MainWindow::saveConfig()
{
	QString basepath(getenv("GADGETRON_HOME"));
	if( basepath==QString("") )
		basepath = QString(".");
	else
		basepath += QString("/config");

	QFileDialog fd;
	fd.setNameFilter(tr("configs (*.xml)"));
	fd.setDirectory(QDir(basepath));
	fd.setFileMode(QFileDialog::AnyFile);
	fd.setAcceptMode(QFileDialog::AcceptSave);

	QStringList fileNames;
	if (fd.exec())
		fileNames = fd.selectedFiles();
	else return;

	QFile file(fileNames[0]);
	if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
		return;
	file.write(plainTextEditConfig->toPlainText().toLatin1());
}

void 
MainWindow::go()
{
	// Let it all happen...
	//
	//

	// Tell Gadgetron which XML gadget configuration script to run.
	if (con_.send_gadgetron_configuration_script(plainTextEditConfig->document()->toPlainText().toStdString()) != 0) {
		ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to send XML configuration to the Gadgetron host")));
		QMessageBox::critical(this, tr("Imageclient"),
									tr("Unable to send XML configuration to the Gadgetron host.\n"),
									QMessageBox::Ok, QMessageBox::Ok);
		return;
	}

	// Pass XML configuration parameters (script defined at the top of the file)
	if (con_.send_gadgetron_parameters(std::string(IMAGE_CLIENT_XML_PARAMETER_SCRIPT)) != 0) {
		ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to send XML parameters to the Gadgetron host")));
		QMessageBox::critical(this, tr("Imageclient"),
											tr("Unable to send XML parameters to the Gadgetron host.\n"),
											QMessageBox::Ok, QMessageBox::Ok);
		return;
	}

	/*
	for (int i = 0; i < num_repetitions_; i++) {



		//We are now ready to send data...
		std::ifstream is;
		is.open (data_file, ios::binary );
		is.seekg (0, ios::end);
		size_t length = is.tellg();
		is.seekg (0, ios::beg);

		while ((length-is.tellg()) > sizeof(GadgetMessageAcquisition)) {
			GadgetContainerMessage<GadgetMessageIdentifier>* m1 =
					new GadgetContainerMessage<GadgetMessageIdentifier>();

			m1->getObjectPtr()->id = GADGET_MESSAGE_ACQUISITION;

			GadgetContainerMessage<GadgetMessageAcquisition>* m2 =
					new GadgetContainerMessage<GadgetMessageAcquisition>();

			//Read acquisition header from disk
			is.read(reinterpret_cast<char*>(m2->getObjectPtr()), sizeof(GadgetMessageAcquisition));

			std::vector<unsigned int> dimensions(2);
			dimensions[0] = m2->getObjectPtr()->samples;
			dimensions[1] = m2->getObjectPtr()->channels;

			GadgetContainerMessage< hoNDArray< std::complex<float> > >* m3 =
					new GadgetContainerMessage< hoNDArray< std::complex< float> > >();

			if (!m3->getObjectPtr()->create(&dimensions)) {
				ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to create storage for NDArray.\n")));
				ACE_DEBUG((LM_ERROR, ACE_TEXT("Requested dimensions were (%d,%d)"), dimensions[0], dimensions[1]));
				return -1;
			}

			//Read data from disk
			is.read(reinterpret_cast<char*>(m3->getObjectPtr()->get_data_ptr()),sizeof(float)*2*m3->getObjectPtr()->get_number_of_elements());

			//Chain the message block together.
			m1->con_t(m2);
			m2->con_t(m3);

			if (con_.putq(m1) == -1) {
				ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to put data package on queue")));
				return -1;
			}
		}

		is.close();

		GadgetContainerMessage<GadgetMessageIdentifier>* m1 =
				new GadgetContainerMessage<GadgetMessageIdentifier>();

		m1->getObjectPtr()->id = GADGET_MESSAGE_CLOSE;

		if (con_.putq(m1) == -1) {
			ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to put CLOSE package on queue")));
			return -1;
		}

		con_.wait();
	}

*/
}

void
MainWindow::set_config(ACE_TCHAR* _filename)
{
	QString basepath(getenv("GADGETRON_HOME"));
	if( basepath==QString("") ) basepath==QString("./");
	else basepath += QString("/config/");

	QString filename = basepath + QString(_filename);
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;

	QByteArray document = file.readAll();
	plainTextEditConfig->setPlainText(document);
}


int
MainWindow::setup_gadgetron_connection()
{
	con_.register_writer(GADGET_MESSAGE_IMAGE_REAL_FLOAT,  new GadgetImageMessageWriter<float>());
	con_.register_reader(GADGET_MESSAGE_IMAGE_REAL_USHORT, new ImageWriter<ACE_UINT16>());
	con_.register_reader(GADGET_MESSAGE_IMAGE_REAL_FLOAT,  new ImageWriter<float>());
	con_.register_reader(GADGET_MESSAGE_IMAGE_CPLX_FLOAT,  new ImageWriter< std::complex<float> >());

	if (con_.open(std::string(hostname_),std::string(port_no_)) != 0) {
		ACE_DEBUG((LM_ERROR, ACE_TEXT("Unable to connect to the Gadgetron host\n")));
		return -1;
	}
	return 0;
}
