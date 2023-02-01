/****************************************************************************

 Copyright (C) 2002-2014 Gilles Debunne. All rights reserved.

 This file is part of the QGLViewer library version 2.6.3.

 http://www.libqglviewer.com - contact@libqglviewer.com

 This file may be used under the terms of the GNU General Public License 
 versions 2.0 or 3.0 as published by the Free Software Foundation and
 appearing in the LICENSE file included in the packaging of this file.
 In addition, as a special exception, Gilles Debunne gives you certain 
 additional rights, described in the file GPL_EXCEPTION in this package.

 libQGLViewer uses dual licensing. Commercial/proprietary software must
 purchase a libQGLViewer Commercial License.

 This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

*****************************************************************************/

#include "qglviewer.h"

#ifndef NO_VECTORIAL_RENDER
#  include "ui_VRenderInterface.h"
#  include "VRender/VRender.h"
#endif

#include "ImageInterface.h"

// Output format list
# include <QImageWriter>

#include <qfileinfo.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qapplication.h>
#include <qmap.h>
#include <qinputdialog.h>
#include <qprogressdialog.h>
#include <qcursor.h>

using namespace std;

////// Static global variables - local to this file //////

static bool  snapshot_format_initialized = false;
// List of available output file formats, formatted for QFileDialog.
static QString formats;
// Converts QFileDialog resulting format to Qt snapshotFormat.
static QMap<QString, QString> Qtformat;
// Converts Qt snapshotFormat to QFileDialog menu string.
static QMap<QString, QString> FDFormatString;
// Converts snapshotFormat to file extension
static QMap<QString, QString> extension;


/*! Sets snapshotFileName(). */
void QGLViewer::setSnapshotFileName(const QString& name)
{
	snapshotFileName_ = QFileInfo(name).absoluteFilePath();
}

#ifndef DOXYGEN
const QString& QGLViewer::snapshotFilename() const
{
	qWarning("snapshotFilename is deprecated. Use snapshotFileName() (uppercase N) instead.");
	return snapshotFileName();
}
#endif


/*! Opens a dialog that displays the different available snapshot formats.

Then calls setSnapshotFormat() with the selected one (unless the user cancels).

Returns \c false if the user presses the Cancel button and \c true otherwise. */
bool QGLViewer::openSnapshotFormatDialog()
{
	bool ok = false;
	QStringList list = formats.split(";;", QString::SkipEmptyParts);
	int current = list.indexOf(FDFormatString[snapshotFormat()]);
	QString format = QInputDialog::getItem(this, "Snapshot format", "Select a snapshot format", list, current, false, &ok);
	if (ok)
		setSnapshotFormat(Qtformat[format]);
	return ok;
}


// Finds all available Qt output formats, so that they can be available in
// saveSnapshot dialog. Initialize snapshotFormat() to the first one.
void QGLViewer::initializeSnapshotFormats()
{
	if (snapshot_format_initialized)
		return;
	QList<QByteArray> list = QImageWriter::supportedImageFormats();
	QStringList formatList;
	for (int i=0; i < list.size(); ++i)
		formatList << QString(list.at(i).toUpper());
	//        qWarning("Available image formats: ");
	//        QStringList::Iterator it = formatList.begin();
	//        while( it != formatList.end() )
	//  	      qWarning((*it++).);  QT4 change this. qWarning no longer accepts QString

#ifndef NO_VECTORIAL_RENDER
	// We add the 3 vectorial formats to the list
	formatList += "EPS";
	formatList += "PS";
	formatList += "XFIG";
#endif

	// Check that the interesting formats are available and add them in "formats"
	// Unused formats: XPM XBM PBM PGM
	QStringList QtText, MenuText, Ext;
	QtText += "PNG";	MenuText += "PNG (*.png)";      Ext += "png";
	QStringList::iterator itText = QtText.begin();
	QStringList::iterator itMenu = MenuText.begin();
	QStringList::iterator itExt  = Ext.begin();

	while (itText != QtText.end())
	{
		//QMessageBox::information(this, "Snapshot ", "Trying format\n"+(*itText));
		if (formatList.contains((*itText)))
		{
			//QMessageBox::information(this, "Snapshot ", "Recognized format\n"+(*itText));
			if (formats.isEmpty())
				setSnapshotFormat(*itText);
			else
				formats += ";;";
			formats += (*itMenu);
			Qtformat[(*itMenu)]  = (*itText);
			FDFormatString[(*itText)]  = (*itMenu);
			extension[(*itText)] = (*itExt);
		}
		// Synchronize parsing
		itText++;
		itMenu++;
		itExt++;
	}
	snapshot_format_initialized = true;
}

// Returns false if the user refused to use the fileName
static bool checkFileName(QString& fileName, QWidget* widget, const QString& snapshotFormat)
{
	if (fileName.isEmpty())
		return false;

	// Check that extension has been provided
	QFileInfo info(fileName);

	if (info.suffix().isEmpty())
	{
		// No extension given. Silently add one
		if (fileName.right(1) != ".")
			fileName += ".";
		fileName += extension[snapshotFormat];
		info.setFile(fileName);
	}
	else if (info.suffix() != extension[snapshotFormat])
	{
		// Extension is not appropriate. Propose a modification
		QString modifiedName = info.absolutePath() + '/' + info.baseName() + "." + extension[snapshotFormat];
		QFileInfo modifInfo(modifiedName);
		int i=(QMessageBox::warning(widget,"Wrong extension",
									info.fileName()+" has a wrong extension.\nSave as "+modifInfo.fileName()+" instead ?",
									QMessageBox::Yes,
									QMessageBox::No,
									QMessageBox::Cancel));
		if (i==QMessageBox::Cancel)
			return false;

		if (i==QMessageBox::Yes)
		{
			fileName = modifiedName;
			info.setFile(fileName);
		}
	}

	return true;
}

#ifndef NO_VECTORIAL_RENDER
// static void drawVectorial(void* param)
void drawVectorial(void* param)
{
	( (QGLViewer*) param )->drawVectorial();
}

#ifndef DOXYGEN
class ProgressDialog
{
public:
	static void showProgressDialog(QGLWidget* parent);
	static void updateProgress(float progress, const QString& stepString);
	static void hideProgressDialog();

private:
	static QProgressDialog* progressDialog;
};

QProgressDialog* ProgressDialog::progressDialog = NULL;

void ProgressDialog::showProgressDialog(QGLWidget* parent)
{
	progressDialog = new QProgressDialog(parent);
	progressDialog->setWindowTitle("Image rendering progress");
	progressDialog->setMinimumSize(300, 40);
	progressDialog->setCancelButton(NULL);
	progressDialog->show();
}

void ProgressDialog::updateProgress(float progress, const QString& stepString)
{
	progressDialog->setValue(int(progress*100));
	QString message(stepString);
	if (message.length() > 33)
		message = message.left(17) + "..." + message.right(12);
	progressDialog->setLabelText(message);
	progressDialog->update();
	qApp->processEvents();
}

void ProgressDialog::hideProgressDialog()
{
	progressDialog->close();
	delete progressDialog;
	progressDialog = NULL;
}

class VRenderInterface: public QDialog, public Ui::VRenderInterface
{
public: VRenderInterface(QWidget *parent) : QDialog(parent) { setupUi(this); }
};

#endif //DOXYGEN

// Pops-up a vectorial output option dialog box and save to fileName
// Returns -1 in case of Cancel, 0 for success and (todo) error code in case of problem.
static int saveVectorialSnapshot(const QString& fileName, QGLWidget* widget, const QString& snapshotFormat)
{
	static VRenderInterface* VRinterface = NULL;

	if (!VRinterface)
		VRinterface = new VRenderInterface(widget);


	// Configure interface according to selected snapshotFormat
	if (snapshotFormat == "XFIG")
	{
		VRinterface->tightenBBox->setEnabled(false);
		VRinterface->colorBackground->setEnabled(false);
	}
	else
	{
		VRinterface->tightenBBox->setEnabled(true);
		VRinterface->colorBackground->setEnabled(true);
	}

	if (VRinterface->exec() == QDialog::Rejected)
		return -1;

	vrender::VRenderParams vparams;
	vparams.setFilename(fileName);

	if (snapshotFormat == "EPS")	vparams.setFormat(vrender::VRenderParams::EPS);
	if (snapshotFormat == "PS")	vparams.setFormat(vrender::VRenderParams::PS);
	if (snapshotFormat == "XFIG")	vparams.setFormat(vrender::VRenderParams::XFIG);

	vparams.setOption(vrender::VRenderParams::CullHiddenFaces, !(VRinterface->includeHidden->isChecked()));
	vparams.setOption(vrender::VRenderParams::OptimizeBackFaceCulling, VRinterface->cullBackFaces->isChecked());
	vparams.setOption(vrender::VRenderParams::RenderBlackAndWhite, VRinterface->blackAndWhite->isChecked());
	vparams.setOption(vrender::VRenderParams::AddBackground, VRinterface->colorBackground->isChecked());
	vparams.setOption(vrender::VRenderParams::TightenBoundingBox, VRinterface->tightenBBox->isChecked());

	switch (VRinterface->sortMethod->currentIndex())
	{
		case 0: vparams.setSortMethod(vrender::VRenderParams::NoSorting); 		break;
		case 1: vparams.setSortMethod(vrender::VRenderParams::BSPSort); 		break;
		case 2: vparams.setSortMethod(vrender::VRenderParams::TopologicalSort); 	break;
		case 3: vparams.setSortMethod(vrender::VRenderParams::AdvancedTopologicalSort);	break;
		default:
			qWarning("VRenderInterface::saveVectorialSnapshot: Unknown SortMethod");
	}

	vparams.setProgressFunction(&ProgressDialog::updateProgress);
	ProgressDialog::showProgressDialog(widget);
	widget->makeCurrent();
	widget->raise();
	vrender::VectorialRender(drawVectorial, (void*) widget, vparams);
	ProgressDialog::hideProgressDialog();
	widget->setCursor(QCursor(Qt::ArrowCursor));

	// Should return vparams.error(), but this is currently not set.
	return 0;
}
#endif // NO_VECTORIAL_RENDER


//class ImageInterface: public QDialog, public Ui::ImageInterface
//{
//public: ImageInterface(QWidget *parent) : QDialog(parent) { setupUi(this); }
//};


// Pops-up an image settings dialog box and save to fileName.
// Returns false in case of problem.
bool QGLViewer::saveImageSnapshot(const QString& fileName)
{
	static ImageInterface* imageInterface = NULL;

	if (!imageInterface)
		imageInterface = new ImageInterface(this);
	else {
		// 	imageInterface->imgScale->setValue(1); // assume the user wants to use the same scale
		int scale = imageInterface->imgScale->value();
		imageInterface->imgWidth->setValue(width() * scale);
		imageInterface->imgHeight->setValue(height() * scale);
		imageInterface->imgQuality->setValue(snapshotQuality());
	}

	if (imageInterface->exec() == QDialog::Rejected)
		return true;

	// Hide closed dialog
	qApp->processEvents();

	setSnapshotQuality(imageInterface->imgQuality->value());

	QColor previousBGColor = backgroundColor();
	if (imageInterface->whiteBackground->isChecked())
		setBackgroundColor(Qt::white);

	QSize finalSize(imageInterface->imgWidth->value(), imageInterface->imgHeight->value());

	qreal oversampling = imageInterface->oversampling->value();
	QSize subSize(int(this->width()/oversampling), int(this->height()/oversampling));

	qreal aspectRatio = width() / static_cast<qreal>(height());
	qreal newAspectRatio = finalSize.width() / static_cast<qreal>(finalSize.height());

	qreal zNear = camera()->zNear();
	qreal zFar = camera()->zFar();

	qreal xMin, yMin;
	bool expand = imageInterface->expandFrustum->isChecked();
	if (camera()->type() == qglviewer::Camera::PERSPECTIVE)
		if ((expand && (newAspectRatio>aspectRatio)) || (!expand && (newAspectRatio<aspectRatio)))
		{
			yMin = zNear * tan(camera()->fieldOfView() / 2.0);
			xMin = newAspectRatio * yMin;
		}
		else
		{
			xMin = zNear * tan(camera()->fieldOfView() / 2.0) * aspectRatio;
			yMin = xMin / newAspectRatio;
		}
	else
	{
		camera()->getOrthoWidthHeight(xMin, yMin);
		if ((expand && (newAspectRatio>aspectRatio)) || (!expand && (newAspectRatio<aspectRatio)))
			xMin = newAspectRatio * yMin;
		else
			yMin = xMin / newAspectRatio;
	}

	QImage image(finalSize.width(), finalSize.height(), QImage::Format_ARGB32);

	if (image.isNull())
	{
		QMessageBox::warning(this, "Image saving error",
							 "Unable to create resulting image",
							 QMessageBox::Ok, QMessageBox::NoButton);
		return false;
	}

	// ProgressDialog disabled since it interfers with the screen grabing mecanism on some platforms. Too bad.
	// ProgressDialog::showProgressDialog(this);

	qreal scaleX = subSize.width() / static_cast<qreal>(finalSize.width());
	qreal scaleY = subSize.height() / static_cast<qreal>(finalSize.height());

	qreal deltaX = 2.0 * xMin * scaleX;
	qreal deltaY = 2.0 * yMin * scaleY;

	int nbX = finalSize.width() / subSize.width();
	int nbY = finalSize.height() / subSize.height();

	// Extra subimage on the right/bottom border(s) if needed
	if (nbX * subSize.width() < finalSize.width())
		nbX++;
	if (nbY * subSize.height() < finalSize.height())
		nbY++;

	makeCurrent();

	// tileRegion_ is used by startScreenCoordinatesSystem to appropriately set the local
	// coordinate system when tiling
	tileRegion_ = new TileRegion();
	qreal tileXMin, tileWidth, tileYMin, tileHeight;
	if ((expand && (newAspectRatio>aspectRatio)) || (!expand && (newAspectRatio<aspectRatio)))
	{
		qreal tileTotalWidth = newAspectRatio * height();
		tileXMin = (width() - tileTotalWidth) / 2.0;
		tileWidth = tileTotalWidth * scaleX;
		tileYMin = 0.0;
		tileHeight = height() * scaleY;
		tileRegion_->textScale = 1.0 / scaleY;
	}
	else
	{
		qreal tileTotalHeight = width() / newAspectRatio;
		tileYMin = (height() - tileTotalHeight) / 2.0;
		tileHeight = tileTotalHeight * scaleY;
		tileXMin = 0.0;
		tileWidth = width() * scaleX;
		tileRegion_->textScale = 1.0 / scaleX;
	}

	int count=0;
	for (int i=0; i<nbX; i++)
		for (int j=0; j<nbY; j++)
		{
			preDraw();

			// Change projection matrix
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			if (camera()->type() == qglviewer::Camera::PERSPECTIVE)
				glFrustum(-xMin + i*deltaX, -xMin + (i+1)*deltaX, yMin - (j+1)*deltaY, yMin - j*deltaY, zNear, zFar);
			else
				glOrtho(-xMin + i*deltaX, -xMin + (i+1)*deltaX, yMin - (j+1)*deltaY, yMin - j*deltaY, zNear, zFar);
			glMatrixMode(GL_MODELVIEW);

			tileRegion_->xMin = tileXMin + i * tileWidth;
			tileRegion_->xMax = tileXMin + (i+1) * tileWidth;
			tileRegion_->yMin = tileYMin + j * tileHeight;
			tileRegion_->yMax = tileYMin + (j+1) * tileHeight;

			draw();
			postDraw();

			// ProgressDialog::hideProgressDialog();
			// qApp->processEvents();

			QImage snapshot = grabFrameBuffer(true);

			// ProgressDialog::showProgressDialog(this);
			// ProgressDialog::updateProgress(count / (qreal)(nbX*nbY),
			// "Generating image ["+QString::number(count)+"/"+QString::number(nbX*nbY)+"]");
			// qApp->processEvents();

			QImage subImage = snapshot.scaled(subSize, Qt::IgnoreAspectRatio, Qt::SmoothTransformation);

			// Copy subImage in image
			for (int ii=0; ii<subSize.width(); ii++)
			{
				int fi = i*subSize.width() + ii;
				if (fi == image.width())
					break;
				for (int jj=0; jj<subSize.height(); jj++)
				{
					int fj = j*subSize.height() + jj;
					if (fj == image.height())
						break;
					image.setPixel(fi, fj, subImage.pixel(ii,jj));
				}
			}
			count++;
		}

	bool saveOK = image.save(fileName, snapshotFormat().toLatin1().constData(), snapshotQuality());

	// ProgressDialog::hideProgressDialog();
	// setCursor(QCursor(Qt::ArrowCursor));

	delete tileRegion_;
	tileRegion_ = NULL;

	if (imageInterface->whiteBackground->isChecked())
		setBackgroundColor(previousBGColor);

	return saveOK;
}


/*! Saves a snapshot of the current image displayed by the widget.

 Options are set using snapshotFormat(), snapshotFileName() and snapshotQuality(). For non vectorial
 image formats, the image size is equal to the current viewer's dimensions (see width() and
 height()). See snapshotFormat() for details on supported formats.

 If \p automatic is \c false (or if snapshotFileName() is empty), a file dialog is opened to ask for
 the file name.

 When \p automatic is \c true, the file name is set to \c NAME-NUMBER, where \c NAME is
 snapshotFileName() and \c NUMBER is snapshotCounter(). The snapshotCounter() is automatically
 incremented after each snapshot saving. This is useful to create videos from your application:
 \code
 void Viewer::init()
 {
   resize(720, 576); // PAL DV format (use 720x480 for NTSC DV)
   connect(this, SIGNAL(drawFinished(bool)), SLOT(saveSnapshot(bool)));
 }
 \endcode
 Then call draw() in a loop (for instance using animate() and/or a camera() KeyFrameInterpolator
 replay) to create your image sequence.

 If you want to create a Quicktime VR panoramic sequence, simply use code like this:
 \code
 void Viewer::createQuicktime()
 {
   const int nbImages = 36;
   for (int i=0; i<nbImages; ++i)
	 {
	   camera()->setOrientation(2.0*M_PI/nbImages, 0.0); // Theta-Phi orientation
	   showEntireScene();
	   update(); // calls draw(), which emits drawFinished(), which calls saveSnapshot()
	 }
 }
 \endcode

 If snapshotCounter() is negative, no number is appended to snapshotFileName() and the
 snapshotCounter() is not incremented. This is useful to force the creation of a file, overwriting
 the previous one.

 When \p overwrite is set to \c false (default), a window asks for confirmation if the file already
 exists. In \p automatic mode, the snapshotCounter() is incremented (if positive) until a
 non-existing file name is found instead. Otherwise the file is overwritten without confirmation.

 The VRender library was written by Cyril Soler (Cyril dot Soler at imag dot fr). If the generated
 PS or EPS file is not properly displayed, remove the anti-aliasing option in your postscript viewer.

 \note In order to correctly grab the frame buffer, the QGLViewer window is raised in front of
 other windows by this method. */
void QGLViewer::saveSnapshot(bool automatic, bool overwrite)
{
	// Ask for file name
	if (snapshotFileName().isEmpty() || !automatic)
	{
		QString fileName;
		QString selectedFormat = FDFormatString[snapshotFormat()];
		fileName = QFileDialog::getSaveFileName(this, "Choose a file name to save under", snapshotFileName(), formats, &selectedFormat,
												overwrite?QFileDialog::DontConfirmOverwrite:QFlags<QFileDialog::Option>(0));
		setSnapshotFormat(Qtformat[selectedFormat]);

		if (checkFileName(fileName, this, snapshotFormat()))
			setSnapshotFileName(fileName);
		else
			return;
	}

	QFileInfo fileInfo(snapshotFileName());

	if ((automatic) && (snapshotCounter() >= 0))
	{
		// In automatic mode, names have a number appended
		const QString baseName = fileInfo.baseName();
		QString count;
		count.sprintf("%.04d", snapshotCounter_++);
		QString suffix;
		suffix = fileInfo.suffix();
		if (suffix.isEmpty())
			suffix = extension[snapshotFormat()];
		fileInfo.setFile(fileInfo.absolutePath()+ '/' + baseName + '-' + count + '.' + suffix);

		if (!overwrite)
			while (fileInfo.exists())
			{
				count.sprintf("%.04d", snapshotCounter_++);
				fileInfo.setFile(fileInfo.absolutePath() + '/' +baseName + '-' + count + '.' + fileInfo.suffix());
			}
	}

	bool saveOK;
#ifndef NO_VECTORIAL_RENDER
	if ( (snapshotFormat() == "EPS") || (snapshotFormat() == "PS") || (snapshotFormat() == "XFIG") )
		// Vectorial snapshot. -1 means cancel, 0 is ok, >0 (should be) an error
		saveOK = (saveVectorialSnapshot(fileInfo.filePath(), this, snapshotFormat()) <= 0);
	else
#endif
		if (automatic)
		{
			QImage snapshot = frameBufferSnapshot();
			saveOK = snapshot.save(fileInfo.filePath(), snapshotFormat().toLatin1().constData(), snapshotQuality());
		}
		else
			saveOK = saveImageSnapshot(fileInfo.filePath());

	if (!saveOK)
		QMessageBox::warning(this, "Snapshot problem", "Unable to save snapshot in\n"+fileInfo.filePath());
}

QImage QGLViewer::frameBufferSnapshot()
{
	// Viewer must be on top of other windows.
	makeCurrent();
	raise();
	// Hack: Qt has problems if the frame buffer is grabbed after QFileDialog is displayed.
	// We grab the frame buffer before, even if it might be not necessary (vectorial rendering).
	// The problem could not be reproduced on a simple example to submit a Qt bug.
	// However, only grabs the backgroundImage in the eponym example. May come from the driver.
	return grabFrameBuffer(true);
}

/*! Same as saveSnapshot(), except that it uses \p fileName instead of snapshotFileName().

 If \p fileName is empty, opens a file dialog to select the name.

 Snapshot settings are set from snapshotFormat() and snapshotQuality().

 Asks for confirmation when \p fileName already exists and \p overwrite is \c false (default).

 \attention If \p fileName is a char* (as is "myFile.jpg"), it may be casted into a \c bool, and the
 other saveSnapshot() method may be used instead. Pass QString("myFile.jpg") as a parameter to
 prevent this. */
void QGLViewer::saveSnapshot(const QString& fileName, bool overwrite)
{
	const QString previousName = snapshotFileName();
	const int previousCounter = snapshotCounter();
	setSnapshotFileName(fileName);
	setSnapshotCounter(-1);
	saveSnapshot(true, overwrite);
	setSnapshotFileName(previousName);
	setSnapshotCounter(previousCounter);
}

/*! Takes a snapshot of the current display and pastes it to the clipboard.

This action is activated by the KeyboardAction::SNAPSHOT_TO_CLIPBOARD enum, binded to \c Ctrl+C by default.
*/
void QGLViewer::snapshotToClipboard()
{
	QClipboard *cb = QApplication::clipboard();
	cb->setImage(frameBufferSnapshot());
}

