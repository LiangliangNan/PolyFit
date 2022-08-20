/*
*  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
*  Copyright (C) 2000-2005 INRIA - Project ALICE
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program; if not, write to the Free Software
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*
*  If you modify this software, you should include a notice giving the
*  name of the person performing the modification, the date of modification,
*  and the reason for such modification.
*
*  Contact: Bruno Levy - levy@loria.fr
*
*     Project ALICE
*     LORIA, INRIA Lorraine,
*     Campus Scientifique, BP 239
*     54506 VANDOEUVRE LES NANCY CEDEX
*     FRANCE
*
*  Note that the GNU General Public License does not permit incorporating
*  the Software into proprietary programs.
*
* As an exception to the GPL, Graphite can be linked with the following (non-GPL) libraries:
*     Qt, SuperLU, WildMagic and CGAL
*/


#include "logger.h"
#include "basic_types.h"
#include "assertions.h"
#include <stdarg.h>

/* 
Disables the warning caused by passing 'this' as an argument while
construction is not finished (in LoggerStream ctor).
As LoggerStreamBuf only stores the pointer for later use, so we can
ignore the fact that 'this' is not completely formed yet.
*/
#pragma warning(disable:4355) 


// Workaround for bugged Wine's built-in MSVCP90.DLL 
#ifdef OGF_OS_WINDOWS_WINE
inline std::ostream& noop(std::ostream& out) { return out ; }
#define OGF_FLUSH noop
#else
#define OGF_FLUSH std::flush
#endif

//_________________________________________________________

int LoggerStreamBuf::sync(){
	std::string str(this->str());
	loggerStream_->notify(str);
	this->str("");
	return 0;
}

//_________________________________________________________

LoggerStream::LoggerStream(Logger* logger)
: std::ostream(new LoggerStreamBuf(this)), logger_(logger) {
}

LoggerStream::~LoggerStream() {
	std::streambuf * buf = rdbuf();
	delete buf;
}

void LoggerStream::notify(std::string& str){
	logger_->notify(this, str);
}


//_________________________________________________________

LoggerClient::~LoggerClient() {
}

//_________________________________________________________

CoutLogger::CoutLogger(){
}

CoutLogger::~CoutLogger(){
	if (Logger::instance())
		Logger::instance()->unregister_client(this);
}

void CoutLogger::out_message(const std::string& value) { std::cout << value << OGF_FLUSH; }
void CoutLogger::warn_message(const std::string& value){ std::cout << value << OGF_FLUSH; }
void CoutLogger::err_message(const std::string& value) { std::cout << value << OGF_FLUSH; }
void CoutLogger::status_message(const std::string& value, int timeout){ }

//_________________________________________________________

FileLogger::FileLogger() : log_file_(nil){
	std::string file_name ;
	if( false ){ // TODO if we can get the file name from environment, set it
		set_file_name( file_name ) ;
	}
}

FileLogger::FileLogger(std::string& file_name) : log_file_(nil){
	set_file_name( file_name ) ;
}

FileLogger::~FileLogger(){
	Logger::instance()->unregister_client(this);
	delete log_file_ ;
	log_file_ = nil ; 
}

void FileLogger::set_file_name(std::string& value){
	log_file_name_ = value ;
	if (log_file_ != nil) {
		delete log_file_ ;
		log_file_ = nil ;
	}
	if(log_file_name_.length() != 0) {
		log_file_ = new std::ofstream(log_file_name_.c_str(), std::ios::app) ;	
	}
}

void FileLogger::out_message(const std::string& value){
	if (log_file_ != nil) {
		*log_file_ << value << std::flush ;
	}
}

void FileLogger::warn_message(const std::string& value){
	if (log_file_ != nil) {
		*log_file_ << value << std::flush ;
	}
}

void FileLogger::err_message(const std::string& value){
	if (log_file_ != nil) {
		*log_file_ << value << std::flush ;
	}
}

void FileLogger::status_message(const std::string& value, int timeout){
}


//_________________________________________________________

Logger* Logger::instance_ = nil ;

void Logger::initialize() {
	instance_ = new Logger() ;
	Logger::out("Logger") << "initialized" << std::endl ;
}

void Logger::terminate() {
	Logger::out("PolyFit") << "terminated... Bye!" << std::endl ;
	ogf_assert(instance_ != nil) ;
	delete instance_ ;
	instance_ = nil ;
}

bool Logger::set_value(FeatureName name, const std::string& value) 
{
	if(name == LOG_FILE_NAME) {
		log_file_name_ = value;
		file_client_ = new FileLogger(log_file_name_);
		register_client(file_client_);
		return true ;
	} 
	else if (name == LOG_REGISTER_FEATURES) {
		std::vector<std::string> features ;
		String::split_string(value,';',features) ;
		log_features_.clear();
		for(unsigned int i=0; i<features.size(); i++) {
			log_features_.insert(features[i]) ;
		}
		log_everything_ = (features.size() == 1 && ((features[0] == "*") || (features[0] == "Everything")));
		return true ;
	} 
	else if (name == LOG_EXCLUDE_FEATURES) {
		std::vector<std::string> features ;
		String::split_string(value,';',features) ;
		log_features_exclude_.clear();
		for(unsigned int i=0; i<features.size(); i++) {
			log_features_exclude_.insert(features[i]) ;
		}
		return true ;
	} 
	else {
		ogf_assert_not_reached;
		return false ;
	}
}

bool Logger::resolve(FeatureName name, std::string& value) const {

	if(name == LOG_FILE_NAME) {
		value = log_file_name_ ;
		return true ;
	}
	else if(name == LOG_REGISTER_FEATURES) {
		value = "" ;
		std::set<std::string>::const_iterator it = log_features_.begin(); 
		for(; it != log_features_.end(); it++) {
			if(value.length() != 0) {
				value += ';' ;
			}
			value += *it ;
		}
		return true ;
	}
	else if(name == LOG_EXCLUDE_FEATURES) {
		value = "" ;
		std::set<std::string>::const_iterator it = log_features_exclude_.begin(); 
		for(; it != log_features_exclude_.end(); it++) {
			if(value.length() != 0) {
				value += ';' ;
			}
			value += *it ;
		}
		return true ;
	} 
	else {
		ogf_assert_not_reached;
		return false ;
	}
}


void Logger::register_client(LoggerClient* c){
	clients.insert(c);
}

void Logger::unregister_client(LoggerClient* c){
	clients.erase(c);
}

bool Logger::is_client(LoggerClient* c){
	return clients.find(c) != clients.end();
}


Logger::Logger() : out_(this), warn_(this), err_(this), status_(this) {
	log_everything_ = false ;

	// add a default client printing stuff to std::cout
	default_client_ = new CoutLogger(); 
	register_client(default_client_ );
	file_client_ = nil ;
}

Logger::~Logger() {
	delete default_client_;
	default_client_ = nil;

	if(file_client_ != nil) {
		delete file_client_;
		file_client_ = nil;
	}
}

LoggerStream& Logger::out(const std::string& feature) {
	ogf_assert(instance_ != nil) ;
	return instance_->out_stream(feature) ;
}

LoggerStream& Logger::err(const std::string& feature) {
	ogf_assert(instance_ != nil) ;
	return instance_->err_stream(feature) ;
}

LoggerStream& Logger::warn(const std::string& feature) {
	ogf_assert(instance_ != nil) ;
	return instance_->warn_stream(feature) ;
}

LoggerStream& Logger::status() {
	ogf_assert(instance_ != nil) ;
	return instance_->status_stream() ;
}

LoggerStream& Logger::out_stream(const std::string& feature) {
	current_feature_ = feature ;
	return out_ ;
}

LoggerStream& Logger::err_stream(const std::string& feature) {
	current_feature_ = feature ;
	return err_ ;
}

LoggerStream& Logger::warn_stream(const std::string& feature) {
	current_feature_ = feature ;
	return warn_ ;
}

LoggerStream& Logger::status_stream() {
	return status_ ;
}


void Logger::notify_out(const std::string& message){
	if( 
		(log_everything_ && log_features_exclude_.find(current_feature_)
		== log_features_exclude_.end() )
		|| (log_features_.find(current_feature_) != log_features_.end())
		) {
			std::set<LoggerClient*>::iterator it = clients.begin();
			if (current_feature_.empty()) {
				for (; it != clients.end(); it++) {
					(*it)->out_message( message );
				}
			} else {
				for (; it != clients.end(); it++) {
					(*it)->out_message( current_feature_ + " " + message );
				}
			}
	} 
}
void Logger::notify_warn(const std::string& message){
	std::set<LoggerClient*>::iterator it = clients.begin();
	if (current_feature_.empty()) {
		for (; it != clients.end(); it++) {
			(*it)->warn_message( "Warning: " + message );
			(*it)->status_message( std::string("Warning: " + message), 1000);
		}
	} 
	else {
		for (; it != clients.end(); it++) {
			(*it)->warn_message( "[" + current_feature_ + "] " + "Warning: " + message );
			(*it)->status_message( std::string("Warning: " + message), 1000);
		}
	}
}
void Logger::notify_err(const std::string& message){
	std::set<LoggerClient*>::iterator it = clients.begin();
	if (current_feature_.empty()) {
		for (; it != clients.end(); it++) {
			(*it)->err_message( "Error: " + message );
			(*it)->status_message( std::string("Error: " + message), 1000);
		}
	}
	else {
		for (; it != clients.end(); it++) {
			(*it)->err_message( "[" + current_feature_ + "] " + "Error: " + message );
			(*it)->status_message( std::string("Error: " + message), 1000);
		}
	}
}
void Logger::notify_status(const std::string& message, int timeout){
	std::set<LoggerClient*>::iterator it;
	for (it = clients.begin(); it != clients.end(); it++){
		(*it)->status_message(message, timeout);
	}
}



void Logger::notify(LoggerStream* s, std::string& message) {
	if(s == &out_) {
		notify_out(message);
	} else if (s == &warn_) {
		notify_warn(message);
	} else if (s == &err_) {
		notify_err(message);
	} else if(s == &status_) {
		notify_status(message, 0);
	} else {
		ogf_assert(false) ;
	}
}


//////////////////////////////////////////////////////////////////////////


void graphite_printf(const char* format, ...) {
	va_list args ;
	char buffer[4096] ;
	va_start(args, format) ;
	vsprintf(buffer, format, args) ;
	va_end(args) ;
	for(char* ptr = buffer; *ptr; ptr++) {
		if(*ptr == '\n') { *ptr = ' '; }
	}
	::Logger::out("") << buffer << std::endl ;
}


void graphite_fprintf(FILE* out, const char* format, ...) {
	va_list args ;
	char buffer[4096] ;
	va_start(args, format) ;
	vsprintf(buffer, format, args) ;
	va_end(args) ;
	for(char* ptr = buffer; *ptr; ptr++) {
		if(*ptr == '\n') { *ptr = ' '; }
	}
	if(out == stdout) {
		::Logger::out("") << buffer << std::endl ;
	} else if(out == stderr) {
		::Logger::err("") << buffer << std::endl ;
	} else {
		fprintf(out, "%s", buffer) ;
	}
}
